from urllib.parse import quote
import argparse
from pathlib import Path
import subprocess as sp

import requests


COMPOUND_DIR_NAME = "pubchem_compounds"
HLGAP_THR = 0.5
VERBOSITY = 3


def rawurlencode(string):
    """
    This function is used to encode the string to be used in the URL.
    """
    return quote(string, safe="-_.~a-zA-Z0-9")


def get_args():
    """
    This function is used to parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Search CIDs from the PubChem Database based on a list of compounds."
    )
    parser.add_argument("compounds", help="Compound or compound list")
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
    return parser.parse_args()


def test_pubchem_server():
    """
    This function is used to test the connection to the PubChem server.
    """
    response = requests.get(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1/cids/TXT", timeout=10
    )
    if response.status_code == 200 and response.text.strip() == "1":
        if VERBOSITY > 1:
            print("PubChem Server is working fine.")
    else:
        raise ConnectionError(
            "No connection could be established. Check if you have access to the internet."
        )


# Function to run the xtb command
def run_xtb_command(xtb_path: Path, calc_dir: Path, args: list) -> tuple:
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


def check_hlgap(output: str) -> bool:
    """
    This function is used to check if the HOMO-LUMO gap of an xtb output is large enough.
    """
    hlgap = None
    for line in output.split("\n"):
        if "HOMO-LUMO GAP" in line:
            hlgap = float(line.split()[3])
    if not hlgap:
        raise ValueError("HOMO-LUMO gap not determined.")

    if hlgap >= HLGAP_THR:
        return True
    else:
        return False


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


def retrieve_3d_sdf(workingdir: Path, cid: str, xtb_path: Path) -> Path:
    """
    This function is used to retrieve the 3D conformer data from the PubChem database.
    """
    response = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/sdf?record_type=3d",
        timeout=10,
    )
    if "PUGREST.NotFound" in response.text:
        try:
            twodthreed(workingdir=workingdir, cid=cid, xtb_path=xtb_path)
        except Exception:
            if VERBOSITY > 1:
                print(f"Error in 2D to 3D structure conversion for CID {cid}.")
            raise Exception
    else:
        with open(workingdir / f"{cid}.sdf", "w", encoding="utf-8") as file:
            file.write(response.text)
        # still conduct single point calculation to get the charge
        xtb_out, xtb_err, exitcode = run_xtb_command(
            xtb_path=xtb_path,
            calc_dir=workingdir,
            args=[f"{cid}.sdf", "--gfn", "2", "--sp", "--ceasefiles"],
        )
        chrg = get_charge_from_xtb_output(xtb_out=xtb_out)
        with open(Path(f"{workingdir}/.CHRG"), "w", encoding="UTF-8") as f:
            f.write(str(chrg) + "\n")
    # return 3d file path
    return workingdir / f"{cid}.sdf"

def get_charge_from_xtb_output(xtb_out: str) -> int:
    """
    This function is used to extract the total charge from the xtb output.
    """
    # load fourth entry of a line with ":: total charge" of xtb.out into a variable
    chrg = 0
    for line in xtb_out.split("\n"):
        if ":: total charge" in line:
            chrg = round(float(line.split()[3]))
            break
    if VERBOSITY > 1:
        print(" " * 3 + f"Total charge: {chrg:6d}")
    return chrg


def twodthreed(workingdir: Path, cid: str, xtb_path: Path):
    if VERBOSITY > 1:
        print(
            f"No 3D Conformer Data found for CID {cid}. Retrieving 2D Conformer Data instead."
        )
    response = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/sdf",
        timeout=10,
    )
    with open(f"{workingdir}/{cid}_2d.sdf", "w", encoding="utf-8") as file:
        file.write(response.text)

    # Convert 2D to 3D
    try:
        xtb_args = [
            f"{cid}_2d.sdf",
            "--gfn",
            "2",
            "--sp",
            "--ceasefiles",
            "--iterations",
            "50",
        ]
        xtb_out, _, exitcode = run_xtb_command(
            xtb_path=xtb_path, calc_dir=workingdir, args=xtb_args
        )
        if exitcode == 0 and "converted geometry written to" in xtb_out:
            if VERBOSITY > 1:
                print("3D conversion successful.")
            Path(f"{workingdir}/gfnff_convert.sdf").rename(f"{workingdir}/{cid}.sdf")
        else:
            raise Exception("3D conversion failed.")
        if not check_hlgap(xtb_out):
            raise ValueError("HOMO-LUMO gap too small!")
        chrg = get_charge_from_xtb_output(xtb_out=xtb_out)
        with open(Path(f"{workingdir}/gfnff_convert.sdf"), "w", encoding="UTF-8") as f:
            f.write(str(chrg) + "\n")
    except Exception as e:
        if VERBOSITY > 1:
            print("xtb structure conversion failed. Original error:")
            print(e)
        raise Exception
    finally:
        files_to_delete = [
            ".sccnotconverged",
            "convert.log",
            "mdrestart",
            "xtb.trj",
            "xtbmdok",
        ]
        for fdel in files_to_delete:
            if Path(f"{workingdir}/{fdel}").exists():
                Path(f"{workingdir}/{fdel}").unlink()


def convert_structure(strucfile: Path, suffix: str) -> Path:
    """
    Conversion of a structure format into another.
    """
    sp.run(
        [
            "mctc-convert",
            strucfile.name,
            strucfile.with_suffix(suffix).name,
            "--normalize",
        ],
        check=True,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        cwd=strucfile.parent,
    )
    return strucfile.with_suffix(suffix)


def opt_structure(xtb_path: Path, strucfile: Path) -> Path:
    xtb_out, xtb_err, exitcode = run_xtb_command(
        xtb_path=xtb_path,
        calc_dir=strucfile.parent,
        args=[strucfile.name, "--opt", "--gfn", "2", "--ceasefiles"],
    )
    # rename the file xtbopt.{suffix} to {strucfile}_opt.{suffix}
    xtb_opt_file = Path(f"{strucfile.parent}/xtbopt{strucfile.suffix}")
    strucfile_opt = Path(f"{strucfile.parent}/{strucfile.stem}_opt{strucfile.suffix}")
    xtb_opt_file.rename(strucfile_opt)
    files_to_delete = [
        ".xtboptok",
    ]
    for fdel in files_to_delete:
        if Path(f"{strucfile.parent}/{fdel}").exists():
            Path(f"{strucfile.parent}/{fdel}").unlink()
    return strucfile_opt

def main():
    """
    This function is used to search the compounds in the PubChem database.
    """
    args = get_args()
    compounds = args.compounds
    input_format = args.input
    output_format = args.output
    skip = args.skip

    xtb_path = "xtb"
    if not skip:
        test_pubchem_server()

    compound_list = []
    if Path(compounds).is_file():
        with open(compounds, "r", encoding="utf-8") as file:
            compound_list = [line.strip() for line in file.readlines()]
    else:
        compound_list = [compounds]

    found_compounds = []
    not_found_compounds = []

    for compound in compound_list:
        result = search_compound(compound, input_format)
        print(result)
        if "PUGREST.NotFound" in result or "PUGREST.BadRequest" in result:
            not_found_compounds.append(compound)
        else:
            found_compounds.append((compound, result))

    if output_format == "sdf":
        for cid, compound in found_compounds:
            cid_dir = Path(f"{COMPOUND_DIR_NAME}/{cid}")
            cid_dir.mkdir(parents=True, exist_ok=True)
            try:
                if VERBOSITY > 2:
                    print(f"Retrieving 3D structure for {compound}.")
                struc = retrieve_3d_sdf(workingdir=cid_dir, cid=cid, xtb_path=xtb_path)
            except Exception as e:
                print(f"Error in retrieving 3D structure for {compound}.")
                print(e)
                not_found_compounds.append(compound)
                continue
            struc = convert_structure(struc, ".xyz")
            try:
                if VERBOSITY > 2:
                    print(f"Optimizing structure for {compound}.")
                struc_opt = opt_structure(xtb_path, struc)
            except Exception as e:
                if VERBOSITY > 1:
                    print(f"Error in optimizing structure for {compound}.")
                    print(e)
                not_found_compounds.append(compound)
                continue
    elif output_format in ["logp", "logP"]:
        with open("pubchem_logP.data", "w", encoding="utf-8") as file:
            for compound, cid in found_compounds:
                response = requests.get(
                    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/XlogP/txt",
                    timeout=10,
                )
                logp = response.text.strip()
                file.write(f"{compound} {cid} {logp}\n")
    elif output_format == "list":
        with open("found.results", "w", encoding="utf-8") as file:
            for compound, cid in found_compounds:
                file.write(f"{compound} {cid}\n")

    if not_found_compounds:
        with open("not_found.compound", "w", encoding="utf-8") as file:
            for compound in not_found_compounds:
                file.write(f"{compound}\n")


if __name__ == "__main__":
    main()
