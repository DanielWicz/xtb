# This file is part of xtb.
#
# Copyright (C) 2019-2021 Sebastian Ehlert
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

# Find test-drive dependency for testing
testdrive_dep = dependency(
  'test-drive',
  version: '>=0.3',
  fallback: ['test-drive', 'testdrive_dep'],
  required: not meson.is_subproject(),
  default_options: [
    'default_library=static',
    'openmp=false',
  ],
)
# If we do not find test-drive and are a subproject, we just skip testing
if not testdrive_dep.found()
  subdir_done()
endif

fs = import('fs')
py = import('python')

tests = [
  'atomlist',
  'coordinationnumber',
  'coulomb',
  'dftd3',
  'dftd4',
  'docking',
  'eeq',
  'gfn0',
  'gfn1',
  'gfn2',
  'gfnff',
  'hessian',
  'palnhess_threads',
  'iff',
  'latticepoint',
  'molecule',
  'oniom',
  'dipro',
  'pbc-tools',
  'peeq',
  'repulsion',
  'symmetry',
  'tblite',
  'thermo',
  'vertical',
  'wsc',
  'cpx',
  'mpi_hessian',
]

test_srcs = files(
  'main.f90',
  'molstock.f90',
  #> PTB as additional test since ending is different
  'test_ptb.F90',
)
foreach t : tests
  fname = 'test_@0@'.format(t.underscorify())
  if fs.is_file(meson.current_source_dir() / (fname + '.F90'))
    test_srcs += files(fname + '.F90')
  else
    test_srcs += files(fname + '.f90')
  endif
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [xtb_dep_static, testdrive_dep],
  fortran_args: [
     '-DWITH_TBLITE=@0@'.format(tblite_dep.found() ? 1 : 0),
     '-DWITH_CPCMX=@0@'.format(cpx_dep.found() ? 1 : 0),
  ],
  link_language: 'fortran',
)

foreach t : tests
  test(t, tester, args: t, suite: 'unit', env: xtbenv)
endforeach
test('ptb', tester, args: 'ptb', suite: 'unit', env: xtbenv)

mpirun_prog = find_program('mpirun', required: false)
if mpirun_prog.found()
  py3 = py.find_installation()
  # mpirun/PMIx in this environment needs AF_INET sockets. Some CI sandboxes
  # forbid them (socket() returns EPERM), which causes mpirun to abort before
  # the test starts. Probe once and skip the mpirun test in that case.
  net_probe = run_command(py3, ['-c', '''
import socket, sys
try:
    socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sys.exit(0)
except OSError:
    sys.exit(1)
'''], check: false)
  mpi_can_launch = net_probe.returncode() == 0

  mpi_env = environment()
  mpi_env.prepend('LD_LIBRARY_PATH', meson.project_build_root())
  mpi_env.prepend('PYTHONPATH', meson.project_source_root() / 'python')
  mpi_env.set('XTBPATH', meson.project_source_root())
  mpi_env.set('PMIX_MCA_ptl', 'tcp')
  mpi_env.set('PMIX_MCA_ptl_tcp_if_include', 'lo')
  mpi_env.set('OMPI_MCA_btl', 'self,tcp')
  mpi_env.set('OMPI_MCA_oob', 'tcp')
  mpi_env.set('OMPI_MCA_rmaps_base_oversubscribe', 'yes')
  if mpi_can_launch
    test('mpi_hessian_mpirun', mpirun_prog, args: ['-np', '3', tester.full_path(), 'mpi_hessian'], suite: 'unit', env: mpi_env)
  else
    message('Skipping mpi_hessian_mpirun: AF_INET sockets not permitted in this environment')
  endif
endif

# some very basic checks to see if the executable reacts at all
test('Argparser print version', xtb_exe, args: '--version', env: xtbenv)
test('Argparser print help', xtb_exe, args: '--help', env: xtbenv)
test('Argparser print license', xtb_exe, args: '--license', env: xtbenv)
test('Argparser no arguments', xtb_exe, should_fail: true, env: xtbenv)
if host_machine.system() != 'windows'
  test('Info', xtb_exe, env: xtbenv, args: [
    'info',
    files('../../assets/inputs/xyz/taxol.xyz'),
    files('../../assets/inputs/coord/caffeine.coord'),
    files('../../assets/inputs/coord/quartz.3d.coord'),
    files('../../assets/inputs/vasp/ammonia.vasp'),
  ])
endif
