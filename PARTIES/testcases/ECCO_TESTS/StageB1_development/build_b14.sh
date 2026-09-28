#!/bin/bash -l
# Roadmap B.1.4 -- build the four binaries the sediment-clip decision needs.
#
#   parties.stageA_b14     Stage-A flags (ECCO code compiled out).  Runs the
#                          recorded Le=100 no-impact gate: must reproduce
#                          lambda = 0.21844955 to every printed digit.
#   parties.stageB_sedctl  2-D Stage-B, SED_CLIP_AUDIT on, restore OFF.
#                          The control: measures the leak directly instead of
#                          inferring it from the meltwater tracer.
#   parties.stageB_sedfix  2-D Stage-B, audit on, restore ON.  Same deck as the
#                          control, one flag apart -- that is the A/B.
#   parties.stageB_sed3d   3-D (x,z periodic, y no-slip), audit on, restore ON,
#                          plus LUBRICATION_NORMAL for the B.2 wall collision.
#
# Boundary.h is restored from a backup at the end, always (trap).
set -euo pipefail

cd /home/x-mjalabert/PARTIES/PARTIES

# The Anvil module init trips over `set -u`, so relax it just for the sourcing.
set +u
source /etc/profile.d/modules.sh
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
set -u

B=src/Include/Boundary.h
BK=$(mktemp /tmp/Boundary.h.b14.XXXXXX)
cp "$B" "$BK"
trap 'cp "$BK" "$B"; rm -f "$BK"; echo "[build_b14] Boundary.h restored"' EXIT

build () {           # build <output-name>
    make clean > /dev/null 2>&1
    make parties > /tmp/b14_build_$1.log 2>&1 || {
        echo "BUILD FAILED: $1 -- tail of log:"; tail -30 /tmp/b14_build_$1.log; exit 1; }
    cp parties "parties.$1"
    printf '  %-22s %10d bytes  md5 %s\n' "$1" \
           "$(stat -c%s parties.$1)" "$(md5sum parties.$1 | cut -c1-8)"
}

on ()  { sed -i "s|^#undef  *$1\b.*|#define $1${2:+ $2}|"  "$B"; }
off () { sed -i "s|^#define  *$1\b.*|#undef $1|"           "$B"; }

echo "[build_b14] 1/4  stageA_b14   (gate: ECCO code compiled out)"
build stageA_b14

echo "[build_b14] 2/4  stageB_sedctl (2-D, audit on, restore OFF)"
on VOF_IBM
on LAG_PARTICLE_RESOLVED
on VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT 200
build stageB_sedctl

echo "[build_b14] 3/5  stageB_sedfix (2-D, audit on, restore ON)"
on VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE
build stageB_sedfix

# ---- ROUND 3: remove the causes instead of accounting for the symptom -------
# Restore OFF from here on: with the CH mask there is nothing to restore, and
# the audit's `net` column is the proof of that (it should collapse to ~0).
off VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE

echo "[build_b14] 4/6  stageB_mask   (2-D, CH mask ON, old shell criterion)"
on VOF_DIFFUSE_SEDIMENT_CH_MASK
on VOF_DIFFUSE_ICE_PENAL_THRESHOLD
build stageB_mask

echo "[build_b14] 5/6  stageB_mask2  (2-D, CH mask ON + non-solid shell norm)"
on VOF_DIFFUSE_SHELL_FRACTION_NONSOLID
build stageB_mask2

echo "[build_b14] 6/6  stageB_sed3d  (3-D, x/z periodic, round-4 no-flux mask)"
off TWOD_CARTESIAN
off PERIODIC_Z_NOSLIP_BOX
on  PERIODIC_NOSLIP_BOX
on  LUBRICATION_NORMAL
build stageB_sed3d

echo "[build_b14] all built."
