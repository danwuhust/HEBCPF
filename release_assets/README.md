# V5.2 release assets

Upload the following six files in the **Assets** area of the HEBCPF v5.2.0
GitHub Release. The filenames retain `v5` because these are the exact verified
V5 artifacts reused by V5.2, not relabeled or regenerated archives:

- `HEBCPF-v5-klurf-SuiteSparse-v7.12.2.mexw64`
- `HEBCPF-v5-klurf-SuiteSparse-v7.12.2.mexw64.sha256`
- `HEBCPF-v5-klurf-source-SuiteSparse-v7.12.2.zip`
- `HEBCPF-v5-klurf-source-SuiteSparse-v7.12.2.zip.sha256`
- `HEBCPF-v5-benchmark-raw.zip`
- `HEBCPF-v5-benchmark-raw.zip.sha256`

## KLU binary and corresponding source

The bundled KLU binary was rebuilt from official SuiteSparse `v7.12.2`, commit
`42151688813c45846a597edcb601435a0e38f3dd`. The source archive contains the
same binary, the HEBCPF wrapper and build helper, exact KLU/BTF/AMD/COLAMD/
SuiteSparse_config source, component licenses, the complete LGPL 2.1 text,
deterministic smoke check, and `BUILD_INFO.md`.

The binary in this asset, `klurf/klurf.mexw64`, and
`HEBCPF_MEX_v5.2/klurf.mexw64` all have SHA-256:

`e84e55c40afd1268a6ccf994df6932b184193fefee33a55bd5523a3a1fee2ff2`

## Scheduler benchmark evidence

`HEBCPF-v5-benchmark-raw.zip` contains the official V5 scheduler dataset from the
designated main job machine. The run used the source-equivalent unreleased V5.5
tree; the authors confirm that V5 differs only by the reproducibly rebuilt KLU
binary. The archive metadata records both KLU hashes, and no result from the
current auxiliary machine is included. This archive is distinct from the V5.2
pure-MATLAB-versus-MEX timing table; raw repeat records for that comparison are
not available and are not claimed to be present in this asset.

## Verification

```powershell
Get-FileHash .\HEBCPF-v5-klurf-SuiteSparse-v7.12.2.mexw64 -Algorithm SHA256
Get-FileHash .\HEBCPF-v5-klurf-source-SuiteSparse-v7.12.2.zip -Algorithm SHA256
Get-FileHash .\HEBCPF-v5-benchmark-raw.zip -Algorithm SHA256
```

Each result must equal the hash stored in the matching `.sha256` file.
