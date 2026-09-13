"""Campaign-selected centered-sieve profile; shared constants with ROOT/C++."""
from dataclasses import dataclass
from pathlib import Path
import argparse
import csv
import math
import re

@dataclass(frozen=True)
class Spectrometer:
    name: str
    prefix: str
    cer: str
    nx: int
    ny: int
    sieve_distance: float
    dx: float
    dy: float
    delta_min: float
    delta_max: float
    xmp: float
    ymp: float
    hb_linear: float
    hb_quadratic: float

    def branch(self, suffix):
        return f"{self.prefix}.{suffix}"

    @property
    def cherenkov_branch(self):
        return self.branch(f"{self.cer}.npeSum")

    def xs(self, index):
        return (index-(self.nx-1)/2)*self.dx

    def ys(self, index):
        return (index-5)*self.dy if self.name=="SHMS" else (index-4)*0.6*2.54


def from_campaign(campaign):
    names = {m.group(1) for part in str(campaign).split("/")
             if (m := re.match(r"^(SHMS|HMS)(?:_|$)", part))}
    if len(names)!=1:
        raise ValueError(f"Expected one unambiguous HMS_<campaign> or SHMS_<campaign> path: {campaign}")
    name=next(iter(names))
    for line in Path(__file__).with_name("spectrometer_profiles.def").read_text().splitlines():
        if not line.startswith("HALLC_PROFILE("):
            continue
        fields=next(csv.reader([line[len("HALLC_PROFILE("):-1]]))
        if fields[0]==name:
            return Spectrometer(*fields[:3], *map(int,fields[3:5]), *map(float,fields[5:]))
    raise ValueError(f"Missing profile: {name}")


def run_metadata(run, filename="DATfiles/list_of_optics_run.dat"):
    """Read target positions and slice edges independently of spectrometer."""
    lines=Path(filename).read_text().splitlines()
    for i,line in enumerate(lines):
        fields=[x.strip() for x in line.split(",")]
        if len(fields)<6 or fields[0]!=str(run):
            continue
        foils=[float(x) for x in lines[i+1].split(",") if x.strip()]
        edges=[float(x) for x in lines[i+2].split(",") if x.strip()]
        if (len(foils)!=int(fields[3]) or len(edges)!=int(fields[5]) or len(edges)<2
            or not all(math.isfinite(x) for x in foils+edges)
            or any(a>=b for a,b in zip(edges,edges[1:]))):
            raise ValueError(f"Invalid foil/delta metadata for {run} in {filename}")
        return dict(angle_deg=float(fields[2]), sieve_flag=int(fields[4]),foils=foils,edges=edges)
    raise ValueError(f"Run {run} not found in {filename}")


def delta_index(edges, low, high):
    matches=[i for i,(a,b) in enumerate(zip(edges,edges[1:]))
             if math.isclose(a,low,abs_tol=1e-8) and math.isclose(b,high,abs_tol=1e-8)]
    if len(matches)!=1:
        raise ValueError(f"Delta interval [{low}, {high}) does not match metadata {edges}")
    return matches[0]


def configure_gmm(args):
    args.spectrometer=spec=from_campaign(args.campaign)
    if args.max_components is None:
        args.max_components=spec.nx if hasattr(args,"yscol") else spec.ny
    for field,default in dict(ysieve_min=-9.0 if spec.name=="SHMS" else -6.8,
                              ysieve_max=9.0 if spec.name=="SHMS" else 6.8,
                              xsieve_min=-14.0 if spec.name=="SHMS" else -12.5,
                              xsieve_max=14.0 if spec.name=="SHMS" else 12.5).items():
        if getattr(args,field) is None:
            setattr(args,field,default)
    return spec


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign")
    parser.add_argument("--run", type=int)
    parser.add_argument("--metadata", default="DATfiles/list_of_optics_run.dat")
    parser.add_argument("--delta-tag")
    args=parser.parse_args()
    spec=from_campaign(args.campaign)
    if args.run is None:
        print(spec.name,spec.name.lower(),spec.nx,spec.ny,f"{spec.delta_min:g}",f"{spec.delta_max:g}",sep="\t")
    else:
        meta=run_metadata(args.run,args.metadata)
        if spec.name=="SHMS" and meta["sieve_flag"]!=1:
            raise ValueError("This version supports centered SHMS sieve only (SieveFlag=1)")
        if args.delta_tag:
            tokens=args.delta_tag.removeprefix("delta_").split("_to_")
            values=[float(t.replace("m","-").replace("p",".")) for t in tokens]
            if len(values)!=2:
                raise ValueError("Expected delta_LOW_to_HIGH tag")
            print(delta_index(meta["edges"],*values))
        else:
            print(" ".join(f"{v:g}" for v in meta["edges"]))

if __name__=="__main__":
    main()
