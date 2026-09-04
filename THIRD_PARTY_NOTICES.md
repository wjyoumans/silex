# Third-Party Notices

## Project license and historical source notice

Silex is distributed under the GNU General Public License, version 3 or, at
your option, any later version (`GPL-3.0-or-later`). The license text is in
[LICENSE](LICENSE).

The initial public Silex source snapshot was assembled from a private
historical `silex-cpp` tree at commit
`7fdce40f30abbb08e024253347efba0d4e5fcc5a`. That identifier is retained as a
provenance record. The private tree has no public remote, so this identifier is
not presented as a publicly reproducible source reference. Reorganizing or
rewriting Git history does not erase upstream authorship, algorithm lineage,
or license obligations. This notice must be retained with redistributed source.

## PARI/GP

Portions of Silex implement or translate algorithms from
[PARI/GP 2.17.3](https://pari.math.u-bordeaux.fr/), including class-group,
unit-group, ideal, relation-search, enumeration, and zeta behavior identified
in the public algorithm/source map.

Copyright (C) 2000-2024 The PARI Group, Bordeaux, with additional authors and
component notices recorded by the PARI/GP distribution.

PARI/GP is free software licensed under the GNU General Public License, version
2 or, at your option, any later version (`GPL-2.0-or-later`). Silex's
`GPL-3.0-or-later` distribution terms preserve the applicable copyleft
obligations for adapted material. Complete PARI/GP author and component notices
remain available from the corresponding PARI/GP source distribution.

PARI/GP is a source-lineage and comparison baseline; Silex does not require a
PARI/GP runtime for its native library.

## Hecke.jl

Selected behavior and algorithm lineage use two fixed public Hecke.jl sources:

- General order, ideal, LLL relation-search, class/unit, compact-element, and
  saturation sources use [Hecke.jl v0.38.6](https://github.com/thofma/Hecke.jl/releases/tag/v0.38.6),
  commit [`74215ba3d34f296e6f709e415e8007d225524287`](https://github.com/thofma/Hecke.jl/commit/74215ba3d34f296e6f709e415e8007d225524287),
  Git tree `5758221d5c6c176b4781dbafe267bb056099d56b`.
- S-unit sources use [Hecke.jl v0.39.19](https://github.com/thofma/Hecke.jl/releases/tag/v0.39.19),
  commit [`122658620f5ac3c8260785c06d9ce7062f037498`](https://github.com/thofma/Hecke.jl/commit/122658620f5ac3c8260785c06d9ce7062f037498),
  Git tree `b538e80d8c6d6eb783b99b0a217c864458cc7b99`.

Hecke.jl is a source and comparison baseline, not a native-library runtime
dependency.

Unless otherwise specified in a file, Hecke is licensed under the BSD 2-Clause
"Simplified" License.

Copyright (c) 2015-2024: Hecke contributors; see the
[Hecke.jl contributors](https://github.com/thofma/Hecke.jl/graphs/contributors).

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Hecke's own notice identifies separately licensed included components; consult
the Hecke.jl distribution when redistributing those components. They are not
bundled in this repository.

## External dependencies

Silex links to FLINT and may optionally use third-party development tools or
backends such as Google Benchmark, fplll, flatter, and Sphinx. Those projects
are not bundled in this repository and remain under the copyright and license
terms supplied by their respective distributions. This file records Silex
source lineage and attribution; it is not a software bill of materials for
every possible toolchain installation.
