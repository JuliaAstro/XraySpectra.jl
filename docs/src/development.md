# Building the documentation

From a local XraySpectra checkout, start Julia with the documentation project:

```bash
cd /Users/james.lee/Documents/WORK/gsoc/XraySpectra
julia --project=docs
```

Set up the local package checkouts:

```julia
import Pkg

Pkg.develop(path=pwd())
Pkg.develop(path="../Spectra.jl")
Pkg.develop(path="../FITSFiles.jl")
Pkg.instantiate()
```

The local directory named `Spectra.jl` currently contains the package named
SpectrumBase.

Build the site from the repository root:

```bash
julia --project=docs docs/make.jl
```

Documenter writes the generated site to `docs/build`. Preview it locally with:

```bash
python3 -m http.server 8000 --directory docs/build
```

Then open [http://localhost:8000](http://localhost:8000). Rebuild and refresh
after changing a documentation source file.
