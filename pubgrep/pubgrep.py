"""
This module contains code for interacting with the power user interface of the PubChem database
and different related functions.
"""

import sys
from urllib.parse import quote
import argparse
from pathlib import Path
import subprocess as sp
import shutil
import multiprocessing as mp
from functools import partial
import random
import time

from tqdm import tqdm
import requests

try:
    from version import __version__  # type: ignore # pylint: disable=import-error
except ImportError:
    __version__ = "0.0.0"


def header() -> str:
    """
    This function is used to print the header of the program.
    """
    headerstr = (
        "---------------------------------------------------------------------\n"
        f"                          PubGrep {__version__}\n"
        "- This Program tries to search CIDs from the Pubchem Database based -\n"
        "- on a list of compounds given as Input. Afterwards it creates sdf  -\n"
        "-   Files for each Compound given in an appropriate subdirectory.   -\n"
        "-     If you are using this program extensively (like, a lot!)      -\n"
        "-   for your Research, please consider citing 10.1039/D3RA01705B    -\n"
        "-                          MS, MM, 2021-2024                        -\n"
        "---------------------------------------------------------------------\n\n"
    )
    return headerstr


# Custom Exceptions for xTB and structure conversions
class XtbFailure(Exception):
    "Raised when the xTB calculation does not provide the expected output."


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
        if verbosity > 3:
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
        # Prepare the data payload
        data = {"inchi": compound}
        search_url = f"{base_url}inchi/{rawurlencode(compound)}/cids/TXT"
    elif input_format == "name":
        search_url = f"{base_url}name/{rawurlencode(compound)}/cids/TXT"
    elif input_format == "cid":
        search_url = f"{base_url}cid/{compound}/property/IUPACname/TXT"
    elif input_format in ["smile", "smiles"]:
        search_url = f"{base_url}smiles/{rawurlencode(compound)}/cids/TXT"
    elif input_format in ["cas", "regid"]:
        search_url = f"{base_url}xref/RegistryID/{rawurlencode(compound)}/cids/TXT"
    else:
        raise ValueError("Invalid input format.")

    if input_format != "inchi":
        response = requests.get(search_url, timeout=10)
    else:
        response = requests.post(search_url, data=data, timeout=10)
    return response.text.strip()


### Class that combines the compound specific data
class Compound:
    """
    This class handles everything related to a compound.
    """

    def __init__(
        self,
        name: str,
        cid: str,
        wdir: Path,
        xtb_path: Path,
        hlgap_thr: float,
        verbosity: int,
    ):
        """
        This function is used to initialize the compound object.
        """
        self.name = name
        self.cid = cid
        self.wdir = wdir
        self.xtb_path: Path = xtb_path
        self.hlgap_thr: float = hlgap_thr
        self.verbosity = verbosity
        self.struc: Path | None = None
        self.chrg: int | None = None
        self.hlgap: float | None = None
        self.logp: float | None = None
        self.numatoms: int | None = None

        # > Create the directory for the compound
        self.wdir.mkdir(parents=True, exist_ok=True)

    def __str__(self):
        comp_string = f"{self.cid}" + f"\t{self.name}"
        if self.struc:
            comp_string += f"\t{self.struc.resolve()}"
        if self.numatoms is not None:
            comp_string += f"\t{self.numatoms}"
        if self.logp is not None:
            comp_string += f"\t{self.logp}"
        if self.chrg is not None:
            comp_string += f"\t{self.chrg}"
        if self.hlgap is not None:
            comp_string += f"\t{self.hlgap:.3f}"

        return comp_string

    def print_csv_header(self):
        """
        This function is used to print the header of the csv file.
        """
        header_string = "CID\tName"
        if self.struc:
            header_string += "\tStructure Path"
        if self.numatoms is not None:
            header_string += "\tNumber of Atoms"
        if self.logp is not None:
            header_string += "\tLogP"
        if self.chrg is not None:
            header_string += "\tTotal Charge"
        if self.hlgap is not None:
            header_string += "\tHOMO-LUMO gap"
        return header_string

    def to_dict(self) -> dict:
        """
        This function is used to convert the compound object to a dictionary.
        """
        comp_dict = {
            "CID": self.cid,
            "Name": self.name,
            "Structure Path": self.struc.resolve() if self.struc else None,
            "Number of Atoms": self.numatoms,
            "LogP": self.logp,
            "Total Charge": self.chrg,
            "HOMO-LUMO gap": self.hlgap,
        }
        return comp_dict

    def retrieve_logp(self):
        """
        This function is used to retrieve the LogP value from the PubChem database.
        """
        response = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{self.cid}/property/XlogP/txt",  # pylint: disable=line-too-long
            timeout=10,
        )
        logp = response.text.strip()
        if "PUGREST.BadRequest" in logp or logp == "":
            raise ValueError(f"LogP not found for {self.cid}.")
        float(logp)
        self.logp = logp

    def retrieve_3d_sdf(self):
        """
        This function is used to retrieve the 3D conformer data from the PubChem database.
        """
        response = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{self.cid}/sdf?record_type=3d",
            timeout=10,
        )
        if "PUGREST.NotFound" in response.text:
            self.twodthreed()
        else:
            with open(self.wdir / f"{self.cid}.sdf", "w", encoding="utf-8") as file:
                file.write(response.text)

        if self.hlgap is None or self.chrg is None:
            xtb_out, _, _ = run_xtb(
                xtb_path=self.xtb_path,
                calc_dir=self.wdir,
                args=[f"{self.cid}.sdf", "--gfn", "2", "--sp", "--ceasefiles"],
            )
            if self.numatoms is None:
                pass
            if self.hlgap is None:
                self.hlgap = get_hlgap_from_xtb_output(xtb_out, self.verbosity)
            if self.chrg is None:
                self.chrg = get_charge_from_xtb_output(xtb_out, self.verbosity)
            with open(Path(f"{self.wdir}/.CHRG"), "w", encoding="UTF-8") as f:
                f.write(str(self.chrg) + "\n")

        if self.hlgap < self.hlgap_thr:
            raise ValueError(
                f"HOMO-LUMO gap too small ({self.hlgap} (is) vs. {self.hlgap_thr} (threshold) eV)"
            )

        self.struc = self.wdir / f"{self.cid}.sdf"

    def twodthreed(self):
        """
        This function is used to convert the 2D structure to 3D structure.
        """
        if self.verbosity > 2:
            print(
                f"No 3D Conformer Data found for CID {self.cid}. "
                + "Retrieving 2D Conformer Data instead."
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
                raise XtbFailure(
                    f"Error in 2D to 3D structure conversion for CID {self.cid}."
                )
            self.hlgap = get_hlgap_from_xtb_output(xtb_out, self.verbosity)
            self.chrg = get_charge_from_xtb_output(
                xtb_out=xtb_out, verbosity=self.verbosity
            )
            with open(Path(f"{self.wdir}/.CHRG"), "w", encoding="UTF-8") as f:
                f.write(str(self.chrg) + "\n")
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
        xtb_out, _, returncode = run_xtb(
            xtb_path=self.xtb_path,
            calc_dir=self.wdir,
            args=[
                self.struc.name,
                "--opt",
                "--gfn",
                "2",
                "--ceasefiles",
            ],
        )

        # delete unnecessary files
        files_to_delete = [
            ".xtboptok",
        ]
        for fdel in files_to_delete:
            if Path(f"{self.wdir}/{fdel}").exists():
                Path(f"{self.wdir}/{fdel}").unlink()
        # define old and new file names
        xtb_opt_file = Path(f"{self.wdir}/xtbopt{self.struc.suffix}")
        strucfile_opt = Path(f"{self.wdir}/{self.struc.stem}_opt{self.struc.suffix}")

        # check if the optimization was successful
        if returncode != 0 or not xtb_opt_file.exists():
            raise XtbFailure("xTB optimization failed.")
        self.hlgap = get_hlgap_from_xtb_output(xtb_out, self.verbosity)
        if self.hlgap < self.hlgap_thr:
            raise ValueError(
                f"HOMO-LUMO gap too small ({self.hlgap} (is) vs. {self.hlgap_thr} (threshold) eV)"
            )

        # rename the optimized structure
        xtb_opt_file.rename(strucfile_opt)
        self.struc = strucfile_opt


def process_compound(
    comp: Compound, optimization: bool = False, verbosity: int = 0
) -> tuple:
    """
    This function is used to process the compound.
    It acts as a wrapper for the multiprocessing pool.
    """
    # > Wait a random time to avoid overloading the server
    time.sleep(random.uniform(0.0, 2.5))

    # > Retrieve the 3D structure
    try:
        if verbosity > 3:
            print(f"Retrieving 3D structure for {comp.cid}.")
        comp.retrieve_3d_sdf()
    except (XtbFailure, ValueError) as e:
        if verbosity > 1:
            print(f"Error in retrieving 3D structure for {comp.cid}. {e}")
        return None, comp

    # > Convert the structure to the desired format
    comp.convert_structure(".xyz")

    # > Get the number of atoms
    with open(comp.struc, "r", encoding="utf-8") as file:  # type: ignore
        xyz = file.read()
    comp.numatoms = get_numatoms_from_xyz(xyz, verbosity)

    # > Optimize the structure
    if optimization:
        try:
            if verbosity > 3:
                print(f"Optimizing structure for {comp.cid}.")
            comp.opt_structure()
        except XtbFailure as e:
            if verbosity > 1:
                print(f"Error in optimizing structure for {comp.cid}. {e}")
            return None, comp

    return comp, None


### Technical functions for running and parsing xtb
def run_xtb(xtb_path: Path, calc_dir: Path, args: list) -> tuple:
    """
    This function is used to run the xtb command.
    """
    singlecore = ["-P", "1"]
    args = singlecore + args
    try:
        xtb_out = sp.run(
            [xtb_path, *args],
            check=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=calc_dir,
            timeout=180,
        )
        stdout = xtb_out.stdout.decode()
        with open(calc_dir / "xtb.out", "w", encoding="utf-8") as file:
            file.write(stdout)
        stderr = xtb_out.stderr.decode()
        with open(calc_dir / "xtb.err", "w", encoding="utf-8") as file:
            file.write(stderr)
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


def get_hlgap_from_xtb_output(output: str, verbosity: int) -> float:
    """
    This function is used to check if the HOMO-LUMO gap of an xtb output is large enough.
    """
    hlgap = None
    for line in output.split("\n"):
        if "HOMO-LUMO GAP" in line:
            hlgap = float(line.split()[3])
            break
    if hlgap is None:
        raise ValueError("HOMO-LUMO gap not determined.")
    if verbosity > 3:
        print(" " * 3 + f"HOMO-LUMO gap: {hlgap:5f}")
    return hlgap


def get_charge_from_xtb_output(xtb_out: str, verbosity: int) -> int:
    """
    This function is used to extract the total charge from the xtb output.
    """
    # load fourth entry of a line with ":: total charge" of xtb.out into a variable
    chrg = None
    for line in xtb_out.split("\n"):
        if ":: total charge" in line:
            chrg = round(float(line.split()[3]))
            break
    if chrg is None:
        raise ValueError("Total charge could not be determined.")
    if verbosity > 3:
        print(" " * 3 + f"Total charge: {chrg:6d}")
    return chrg


def get_numatoms_from_xyz(xyz: str, verbosity: int) -> int:
    """
    This function is used to extract the total charge from the xtb output.
    """
    nat = None
    # get number of atoms from first line of XYZ file
    nat = int(xyz.split("\n")[0])
    if verbosity > 3:
        print(" " * 3 + f"Total number of atoms: {nat}")
    return nat


### Argument parsing functions
def get_args() -> argparse.Namespace:
    """
    This function is used to parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description=header()
        + "\nThis program searches CIDs from the PubChem Database based on a list of compounds.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
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
    try:
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
        return 0
    except (ValueError, XtbFailure) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


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
) -> list[Compound]:
    """
    This function is used to search the compounds in the PubChem database.
    """

    if skip:
        test_pubchem_server(verbosity)

    # if len of compounds is 1, take it as a single compound
    if isinstance(compounds, list):
        if len(compounds) == 1:
            compounds = compounds[0]

    compound_list = []
    if isinstance(compounds, list):
        for compound in compounds:
            compound_list.append(compound)
    elif isinstance(compounds, Path) or Path(compounds).is_file():  # type: ignore
        with open(compounds, "r", encoding="utf-8") as file:  # type: ignore
            compound_list = [line.strip() for line in file.readlines()]
    elif isinstance(compounds, str):
        compound_list.append(compounds)
    else:
        raise ValueError("Invalid input for compounds.")

    if len(compound_list) == 0:
        raise ValueError("No compounds found in the input file.")

    found_compounds: list[Compound] = []
    not_found_compounds: list[str] = []

    for compound in compound_list:
        result = search_compound(compound, input_format)
        if input_format != "cid":
            cid = result
            name = compound
        else:
            cid = compound
            name = result
        if verbosity > 3:
            print(f"Found CID {cid} for {name}.")
        if "PUGREST.NotFound" in result or "PUGREST.BadRequest" in result:
            not_found_compounds.append(compound)
        else:
            comp = Compound(
                name=name,
                cid=cid,
                wdir=Path(basedir / f"{cid}"),
                xtb_path=xtb_path,
                hlgap_thr=hlgap_thr,
                verbosity=verbosity,
            )
            found_compounds.append(comp)
            if verbosity > 1:
                print(comp)

    failed_compounds = []
    successful_compounds = []

    if output_format == "list":
        with open("found_compounds.csv", "w", encoding="utf-8") as file:
            for comp in found_compounds:
                print(comp, file=file)
    elif output_format == "sdf":
        wrap_process_compound = partial(
            process_compound, optimization=optimization, verbosity=verbosity
        )
        if verbosity > 3:
            print(f"Using {mp.cpu_count()} cores for processing.")
        with mp.Pool(processes=mp.cpu_count()) as pool:
            for result in tqdm(
                pool.imap(wrap_process_compound, found_compounds),
                total=len(found_compounds),
            ):
                success, failure = result
                if success:
                    successful_compounds.append(success)
                if failure:
                    failed_compounds.append(failure)

    elif output_format in ["logp", "logP"]:
        for comp in found_compounds:
            try:
                if verbosity > 3:
                    print(f"Retrieving LogP for {comp.cid}.")
                comp.retrieve_logp()
            except ValueError as e:
                if verbosity > 1:
                    print(f"Error in retrieving LogP for {comp.cid}. {e}")
                failed_compounds.append(comp)
                continue
            successful_compounds.append(comp)
    else:
        raise ValueError("Invalid output format.")

    if successful_compounds:
        successful_compounds.sort(key=lambda x: int(x.cid))
        with open("compounds.csv", "w", encoding="utf-8") as file:
            print(successful_compounds[0].print_csv_header(), file=file)
            for comp in successful_compounds:
                print(comp, file=file)
    else:
        raise ValueError("No compounds could be processed.")

    if not_found_compounds:
        with open("not_found.compound", "w", encoding="utf-8") as file:
            for compound in not_found_compounds:
                file.write(f"{compound}\n")

    return successful_compounds


if __name__ == "__main__":
    raise SystemExit(cli())
