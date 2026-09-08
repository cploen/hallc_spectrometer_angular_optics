"""Synthetic ROOT checks only; these do not validate the replay data or optics."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


REPO = Path(__file__).resolve().parents[1]
MACRO = REPO / "diagnostics/validation/replay/probe_replay_inputs.C"


def root(*expressions):
    command = ["root", "-l", "-b", "-q"]
    for expression in expressions:
        command.extend(["-e", expression])
    return subprocess.run(command, text=True, capture_output=True, timeout=90)


def probe(campaign, path):
    return root(f'#include {json.dumps(str(MACRO))}',
                f'gSystem->Exit(probe_replay_inputs({json.dumps(campaign)},'
                f'{json.dumps(str(path))},50));')


def main():
    with tempfile.TemporaryDirectory(prefix='hallc-probe-test-') as temp:
        folder = Path(temp)
        generator = folder / "fixture.C"
        generator.write_text(r'''
#include <TFile.h>
#include <TTree.h>
#include <map>
#include <string>
#include <limits>
void fixture(const char* path, int mode) {
  TFile file(path, "RECREATE");
  TTree tree("T", "Synthetic replay; not physics data");
  std::map<std::string, double> values;
  for (const auto* arm : {"H", "P"}) {
    std::string prefix = std::string(arm) + ".";
    for (const auto* suffix : {"cal.etottracknorm", "gtr.dp", "gtr.x", "gtr.y",
          "gtr.th", "gtr.ph", "dc.x_fp", "dc.xp_fp", "dc.y_fp", "dc.yp_fp",
          "react.x", "react.y", "react.z", "extcor.xsieve", "extcor.ysieve",
          "rb.raster.fr_xbpm_tar", "rb.raster.fr_ybpm_tar"}) values[prefix+suffix]=0;
    if (!(mode==1 && prefix=="P."))
      values[prefix+(prefix=="P." ? "ngcer.npeSum" : "cer.npeSum")]=0;
  }
  // Supplied replay schema: react branches present, P BPM absent.
  values.erase("P.rb.raster.fr_xbpm_tar");
  values.erase("P.rb.raster.fr_ybpm_tar");
  for (auto& item : values)
    tree.Branch(item.first.c_str(), &item.second, (item.first+"/D").c_str());
  double dp[] = {-12,-5,0,8,15,0};
  double npe[] = {7,7,3,7,7,7};
  double cal[] = {.9,.7,.9,.9,.9,.8};
  for (int i=0; i<6; ++i) {
    values["P.gtr.dp"] = dp[i];
    values["P.cal.etottracknorm"] = cal[i];
    if (mode!=1) values["P.ngcer.npeSum"] = npe[i];
    if (mode==2 && i==0) values["P.react.x"] = std::numeric_limits<double>::quiet_NaN();
    else values["P.react.x"] = 0;
    if (mode!=3) tree.Fill();
  }
  tree.Write();
}
''')
        cases = [folder / f"case_{i}.root" for i in range(4)]
        expressions = [f'#include {json.dumps(str(generator))}']
        expressions += [f'fixture({json.dumps(str(path))},{i});'
                        for i, path in enumerate(cases)]
        result = root(*expressions)
        assert result.returncode == 0, result.stdout + result.stderr
        checksums = [hashlib.sha256(path.read_bytes()).hexdigest() for path in cases]
        result = probe("/work/SHMS_8p5695GeV/06a_fit_ntuple/root/", cases[0])
        assert result.returncode == 0, result.stdout + result.stderr
        expected = {"ridge_npe>2_-10<dp<10": 4,
                    "npe>6_cal>0.65_no_dp_window": 5,
                    "npe>6_cal>0.65_-10<dp<10": 3,
                    "npe>6_cal>0.65_-10<dp<22": 4,
                    "npe>6_cal>0.8_no_dp_window": 3,
                    "npe>6_cal>0.8_-10<dp<10": 1,
                    "npe>6_cal>0.8_-10<dp<22": 2,
                    "npe>6_cal>0.8_-15<dp<24": 3}
        counts = {line.split("\t")[1]: int(line.split("\t")[2])
                  for line in result.stdout.splitlines() if line.startswith("CUT\t")}
        assert counts == expected, counts
        assert "BRANCH\tH." not in result.stdout
        result = probe("HMS_test", cases[0])
        assert result.returncode == 0, result.stdout + result.stderr
        assert all(line.endswith("\t0") for line in result.stdout.splitlines()
                   if line.startswith("CUT\t")), result.stdout
        for campaign in ("unknown", "HMSfoo", "SHMSfoo", "HMS_outer/SHMS_inner"):
            assert probe(campaign, cases[0]).returncode == 2, campaign
        for path in cases[1:]:
            result = probe("SHMS_test", path)
            assert result.returncode == 3, result.stdout + result.stderr
        assert probe("SHMS_test", folder / "absent.root").returncode == 2
        # End-to-end wrapper, including shell/C++ escaping in replay paths.
        replay_dir = folder / 'replays with "quotes" and spaces'
        replay_dir.mkdir()
        for run in (3283, 3284, 3285, 3286):
            shutil.copy2(cases[0], replay_dir / f"deut_replay_prod_{run}_-1.root")
        env = dict(os.environ, REPLAY_INPUT_DIR=str(replay_dir), PROBE_MAX_EVENTS="6")
        result = subprocess.run(["bash", str(REPO / "diagnostics/run_replay_probe.sh")],
                                env=env, text=True, capture_output=True, timeout=120)
        assert result.returncode == 0, result.stdout + result.stderr
        assert result.stdout.count("STATUS\tSCHEMA_SAMPLE_OK") == 4
        assert checksums == [hashlib.sha256(path.read_bytes()).hexdigest() for path in cases]
        print("PASS: campaign routing, cut counts, missing/nonfinite/empty inputs, wrapper, unchanged input files")


if __name__ == "__main__":
    main()
