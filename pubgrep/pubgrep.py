from urllib.parse import quote
import argparse
from pathlib import Path
import subprocess as sp
import shutil

import requests

try:
    from ._version import version as __version__ # type: ignore
except ImportError:
    __version__ = "0.0.0"

### URL request functions


def rawurlencode(string):
    """
    This function is used to encode the string to be used in the URL.
    """
    return quote(string, safe="-_.~a-zA-Z0-9")


def test_pubchem_server(verbosity):
    """
    This function is used to test the connection to the PubChem server.
    """
    response = requests.get(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1/cids/TXT", timeout=10
    )
    if response.status_code == 200 and response.text.strip() == "1":
        if verbosity > 1:
            print("PubChem Server is working fine.")
    else:
        raise ConnectionError(
            "No connection could be established. Check if you have access to the internet."
        )


def search_compound(compound, input_format):
    """
    This function is used to search the compound in the PubChem database.
    """
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
    search_url = ""
    if input_format == "inchi":
        search_url = f"{base_url}inchi/{rawurlencode(compound)}/cids/TXT"
    elif input_format == "name":
        search_url = f"{base_url}name/{rawurlencode(compound)}/cids/TXT"
    elif input_format == "cid":
        search_url = f"{base_url}cid/{compound}/property/IUPACname/TXT"
    elif input_format in ["smile", "smiles"]:
        search_url = f"{base_url}smiles/{rawurlencode(compound)}/cids/TXT"
    elif input_format in ["cas", "regid"]:
        search_url = f"{base_url}xref/RegistryID/{rawurlencode(compound)}/cids/TXT"

    response = requests.get(search_url, timeout=10)
    return response.text.strip()


### Class that combines the compound specific data


class Compound:
    """
    This class handles everything related to a compound.
    """

    HLGAP_THR: float = 0.5

    def __init__(self, name: str, cid: str, wdir: Path, xtb_path: Path, verbosity: int):
        """
        This function is used to initialize the compound object.
        """
        self.name = name
        self.cid = cid
        self.wdir = wdir
        self.xtb_path: Path = xtb_path
        self.verbosity = verbosity
        self.struc: Path | None = None
        self.chrg: int | None = None
        self.hlgap: float | None = None

    def __str__(self):
        return f"{self.cid}\t{self.name}{"\t"+str(self.struc.resolve()) if self.struc else ''}"

    def retrieve_3d_sdf(self):
        """
        This function is used to retrieve the 3D conformer data from the PubChem database.
        """
        response = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{self.cid}/sdf?record_type=3d",
            timeout=10,
        )
        if "PUGREST.NotFound" in response.text:
            try:
                self.twodthreed()
            except Exception as exc:
                if self.verbosity > 1:
                    print(f"Error in 2D to 3D structure conversion for CID {self.cid}.")
                raise Exception from exc
        else:
            with open(self.wdir / f"{self.cid}.sdf", "w", encoding="utf-8") as file:
                file.write(response.text)

        if not self.hlgap or not self.chrg:
            xtb_out, _, _ = run_xtb(
                xtb_path=self.xtb_path,
                calc_dir=self.wdir,
                args=[f"{self.cid}.sdf", "--gfn", "2", "--sp", "--ceasefiles"],
            )
            self.hlgap = get_hlgap_from_xtb_output(xtb_out)
            self.chrg = get_charge_from_xtb_output(xtb_out, self.verbosity)
            with open(Path(f"{self.wdir}/.CHRG"), "w", encoding="UTF-8") as f:
                f.write(str(self.chrg) + "\n")

        if self.hlgap < self.HLGAP_THR:
            raise ValueError("HOMO-LUMO gap too small!")

        self.struc = self.wdir / f"{self.cid}.sdf"

    def twodthreed(self):
        """
        This function is used to convert the 2D structure to 3D structure.
        """
        if self.verbosity > 1:
            print(
                f"No 3D Conformer Data found for CID {self.cid}. Retrieving 2D Conformer Data instead."
            )
        response = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{self.cid}/sdf",
            timeout=10,
        )
        with open(f"{self.wdir}/{self.cid}_2d.sdf", "w", encoding="utf-8") as file:
            file.write(response.text)

        # Convert 2D to 3D
        try:
            xtb_args = [
                f"{self.cid}_2d.sdf",
                "--gfn",
                "2",
                "--sp",
                "--ceasefiles",
                "--iterations",
                "50",
            ]
            xtb_out, _, exitcode = run_xtb(
                xtb_path=self.xtb_path, calc_dir=self.wdir, args=xtb_args
            )
            if exitcode == 0 and "converted geometry written to" in xtb_out:
                if self.verbosity > 1:
                    print("3D conversion successful.")
                Path(f"{self.wdir}/gfnff_convert.sdf").rename(
                    f"{self.wdir}/{self.cid}.sdf"
                )
            else:
                raise Exception("3D conversion failed.")
            self.hlgap = get_hlgap_from_xtb_output(xtb_out)
            self.chrg = get_charge_from_xtb_output(
                xtb_out=xtb_out, verbosity=self.verbosity
            )
            with open(
                Path(f"{self.wdir}/gfnff_convert.sdf"), "w", encoding="UTF-8"
            ) as f:
                f.write(str(self.chrg) + "\n")
        except Exception as e:
            if self.verbosity > 1:
                print(f"xtb structure conversion failed. Original error:\n{e}")
            raise Exception from e
        finally:
            files_to_delete = [
                ".sccnotconverged",
                "convert.log",
                "mdrestart",
                "xtb.trj",
                "xtbmdok",
            ]
            for fdel in files_to_delete:
                if Path(f"{self.wdir}/{fdel}").exists():
                    Path(f"{self.wdir}/{fdel}").unlink()

    def convert_structure(self, suffix: str):
        """
        Conversion of a structure format into another.
        """
        if self.struc:
            sp.run(
                [
                    "mctc-convert",
                    self.struc.name,
                    self.struc.with_suffix(suffix).name,
                    "--normalize",
                ],
                check=True,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                cwd=self.struc.parent,
            )
            self.struc = self.struc.with_suffix(suffix)

    def opt_structure(self):
        """
        This function is used to optimize the structure of the compound.
        """
        _, _, _ = run_xtb(
            xtb_path=self.xtb_path,
            calc_dir=self.wdir,
            args=[self.struc.name, "--opt", "--gfn", "2", "--ceasefiles"],
        )
        xtb_opt_file = Path(f"{self.wdir}/xtbopt{self.struc.suffix}")
        strucfile_opt = Path(f"{self.wdir}/{self.struc.stem}_opt{self.struc.suffix}")
        xtb_opt_file.rename(strucfile_opt)
        files_to_delete = [
            ".xtboptok",
        ]
        for fdel in files_to_delete:
            if Path(f"{self.wdir}/{fdel}").exists():
                Path(f"{self.wdir}/{fdel}").unlink()
        self.struc = strucfile_opt


### Technical functions for running and parsing xtb


def run_xtb(xtb_path: Path, calc_dir: Path, args: list) -> tuple:
    """
    This function is used to run the xtb command.
    """
    try:
        xtb_out = sp.run(
            [xtb_path, *args],
            check=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=calc_dir,
        )
        stdout = xtb_out.stdout.decode()
        stderr = xtb_out.stderr.decode()
        returncode = xtb_out.returncode
    except sp.CalledProcessError as e:
        with open(calc_dir / "xtb.out", "w", encoding="utf-8") as file:
            file.write(e.stdout.decode())
        with open(calc_dir / "xtb.err", "w", encoding="utf-8") as file:
            print(
                f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.",
                file=file,
            )
            file.write(e.stderr.decode())
        stdout = e.stdout.decode()
        stderr = e.stderr.decode()
        returncode = e.returncode

    return stdout, stderr, returncode


def get_hlgap_from_xtb_output(output: str) -> float:
    """
    This function is used to check if the HOMO-LUMO gap of an xtb output is large enough.
    """
    hlgap = None
    for line in output.split("\n"):
        if "HOMO-LUMO GAP" in line:
            hlgap = float(line.split()[3])
    if not hlgap:
        raise ValueError("HOMO-LUMO gap not determined.")
    return hlgap


def get_charge_from_xtb_output(xtb_out: str, verbosity: int) -> int:
    """
    This function is used to extract the total charge from the xtb output.
    """
    # load fourth entry of a line with ":: total charge" of xtb.out into a variable
    chrg = 0
    for line in xtb_out.split("\n"):
        if ":: total charge" in line:
            chrg = round(float(line.split()[3]))
            break
    if verbosity > 1:
        print(" " * 3 + f"Total charge: {chrg:6d}")
    return chrg


### Argument parsing functions


def get_args() -> argparse.Namespace:
    """
    This function is used to parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Search CIDs from the PubChem Database based on a list of compounds."
    )
    # Required arguments (no flag but positional arguments, doesn't matter how many)
    parser.add_argument(
        "compounds",
        help="List of compounds to search for.",
        nargs="+",
        type=str,
        default=None,
        action="store",
    )
    parser.add_argument(
        "--input",
        default="name",
        choices=["name", "cid", "smiles", "cas", "inchi"],
        help="Input format",
    )
    parser.add_argument(
        "--output", default="sdf", choices=["sdf", "logP", "list"], help="Output format"
    )
    parser.add_argument(
        "--skip", action="store_true", help="Skip PubChem server testing"
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        choices=[0, 1, 2, 3],
        help="Verbosity level",
        required=False,
    )
    parser.add_argument(
        "--opt",
        action="store_true",
        help="Optimize the structure of the compound",
        required=False,
        default=False,
    )
    parser.add_argument(
        "--basedir",
        type=Path,
        default="pubchem_compounds",
        help="Base directory for the compound data",
        required=False,
    )
    parser.add_argument(
        "--hlgap_thr",
        type=float,
        default=0.5,
        help="Threshold for the HOMO-LUMO gap",
        required=False,
    )
    return parser.parse_args()


def cli():
    """
    This function is used to parse the command line arguments.
    """
    args = get_args()
    verbosity = args.verbosity
    compounds = args.compounds
    input_format = args.input
    output_format = args.output
    optimization = args.opt
    skip = args.skip
    basedir = args.basedir
    hlgap_thr = args.hlgap_thr

    xtb_path = Path(shutil.which("xtb")).resolve()

    pubgrep(
        xtb_path=xtb_path,
        basedir=basedir,
        compounds=compounds,
        input_format=input_format,
        output_format=output_format,
        optimization=optimization,
        hlgap_thr=hlgap_thr,
        skip=skip,
        verbosity=verbosity,
    )


def pubgrep(
    xtb_path: Path,
    basedir: Path,
    compounds: list[str] | Path | str,
    input_format: str,
    output_format: str,
    optimization: bool,
    hlgap_thr: float,
    skip: bool,
    verbosity: int,
):
    """
    This function is used to search the compounds in the PubChem database.
    """

    if output_format not in ["sdf", "logP", "list"]:
        raise ValueError("Invalid output format.")

    if not skip:
        test_pubchem_server(verbosity)

    compound_list = []
    if isinstance(compounds, list):
        for compound in compounds:
            compound_list.append(compound)
    if isinstance(compounds, Path):
        if Path(compounds).is_file():
            with open(compounds, "r", encoding="utf-8") as file:
                compound_list = [line.strip() for line in file.readlines()]
    if isinstance(compounds, str):
        compound_list.append(compounds)

    if compound_list == []:
        raise ValueError("No compounds found in the input file.")

    found_compounds: list[Compound] = []
    not_found_compounds: list[str] = []
    Compound.HLGAP_THR = hlgap_thr

    for compound in compound_list:
        result = search_compound(compound, input_format)
        if "PUGREST.NotFound" in result or "PUGREST.BadRequest" in result:
            not_found_compounds.append(compound)
        else:
            comp = Compound(
                name=result,
                cid=compound,
                wdir=Path(basedir / f"{compound}"),
                xtb_path=xtb_path,
                verbosity=verbosity,
            )
            found_compounds.append(comp)
            if verbosity > 0:
                print(comp)

    if output_format == "list":
        with open("found_compounds.csv", "w", encoding="utf-8") as file:
            for comp in found_compounds:
                print(comp, file=file)
    elif output_format == "sdf":
        failed_compounds = []
        successful_compounds = []
        for comp in found_compounds:
            comp.wdir.mkdir(parents=True, exist_ok=True)
            try:
                if verbosity > 2:
                    print(f"Retrieving 3D structure for {comp}.")
                comp.retrieve_3d_sdf()
            except Exception as e:
                if verbosity > 1:
                    print(f"Error in retrieving 3D structure for {comp}.\n{e}")
                failed_compounds.append(comp)
                continue
            comp.convert_structure(".xyz")
            if optimization:
                try:
                    if verbosity > 2:
                        print(f"Optimizing structure for {comp}.")
                    comp.opt_structure()
                except Exception as e:
                    if verbosity > 1:
                        print(f"Error in optimizing structure for {comp}.\n{e}")
                    failed_compounds.append(comp)
                    continue
            successful_compounds.append(comp)

        if successful_compounds:
            with open("compounds.csv", "w", encoding="utf-8") as file:
                for comp in successful_compounds:
                    print(comp, file=file)
    elif output_format in ["logp", "logP"]:
        raise NotImplementedError("LogP calculation is not implemented yet.")
        # with open("pubchem_logP.data", "w", encoding="utf-8") as file:
        #     for compound, cid in found_compounds:
        #         response = requests.get(
        #             f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/XlogP/txt",
        #             timeout=10,
        #         )
        #         logp = response.text.strip()
        #         file.write(f"{compound} {cid} {logp}\n")

    if not_found_compounds:
        with open("not_found.compound", "w", encoding="utf-8") as file:
            for compound in not_found_compounds:
                file.write(f"{compound}\n")


if __name__ == "__main__":
    cli()