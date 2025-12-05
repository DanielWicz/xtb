#!/usr/bin/env python3
"""
Helper wrapper to run the mpi_hessian test under mpirun.

Exits with status 77 (Meson skip) if AF_INET sockets are not permitted
or if mpirun fails with a PMIx listener error, so CI environments without
network sockets don't fail the suite.
"""

import os
import socket
import subprocess
import sys

MPIRUN = sys.argv[1]
TESTER = sys.argv[2]


def can_use_inet() -> bool:
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.close()
        return True
    except OSError:
        return False


def main() -> int:
    if not can_use_inet():
        print("Skipping mpi_hessian_mpirun: AF_INET sockets not permitted")
        return 77

    cmd = [MPIRUN, "-np", "3", TESTER, "mpi_hessian"]
    try:
        completed = subprocess.run(
            cmd, check=False, stdout=sys.stdout, stderr=sys.stderr
        )
    except FileNotFoundError:
        print("Skipping mpi_hessian_mpirun: mpirun not found")
        return 77

    rc = completed.returncode
    # OpenMPI returns 213 (signal 85) on PMIx listener failure in restricted envs.
    if rc in (213, 85):
        print("Skipping mpi_hessian_mpirun: mpirun failed to start PMIx listener")
        return 77
    return rc


if __name__ == "__main__":
    sys.exit(main())
