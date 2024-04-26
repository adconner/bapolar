let
  pkgs = import (fetchTarball "https://github.com/NixOS/nixpkgs/tarball/nixos-23.11") { config = {}; overlays = []; };
  
  papilo = pkgs.stdenv.mkDerivation {
    pname = "papilo";
    version = "2.2.0";
    src = pkgs.fetchFromGitHub {
      owner = "scipopt";
      repo = "papilo";
      rev = "v2.2.0";
      hash = "sha256-X6xr7nhTj5q8QJHn4AtUZSTyVusUDv5X4Dgv0bLf0kE=";
    };
    nativeBuildInputs = with pkgs; [
      cmake
      gfortran
    ];
    buildInputs = with pkgs; [
      boost
      tbb
      gmp
      blas
      lapack
    ];
  };
  soplex = pkgs.stdenv.mkDerivation {
  # currently seems to not link papilo, why?
    pname = "soplex";
    version = "7.0.0";
    src = pkgs.fetchFromGitHub {
      owner = "scipopt";
      repo = "soplex";
      rev = "release-700";
      hash = "sha256-biy69IjqVncvdPhYP83odiDb+1AO0NOMxe/dDEetuJE=";
    };
    # cmakeFlags = [ "-DPAPILO=on" ];
    nativeBuildInputs = with pkgs; [
      cmake
    ];
    buildInputs = with pkgs; [
      gmp
      boost
      tbb
      blas
      lapack
      #papilo
      zlib
      mpfr
      ];
  };
  scip = pkgs.stdenv.mkDerivation {
    pname = "scip";
    version = "9.0.0";
    src = pkgs.fetchFromGitHub {
      owner = "scipopt";
      repo = "scip";
      rev = "v900";
      hash = "sha256-V3NaYoH+4GIv4MWJwgmqIJTfEdsKnPxaYLl7qznejp0=";
    };

    cmakeFlags = [ "-DAUTOBUILD=on" ];
    nativeBuildInputs = with pkgs; [
      cmake
      wget
      unzip
      git
      pkg-config
      flex
      bison
      gfortran
      file
    ];
    buildInputs = with pkgs; [
      gmp
      gnum4
      xz
      zlib
      readline
      boost
      tbb
      blas
      lapack
      gsl
      cliquer
      metis
      #hmetis
      ipopt
      mpfr
      soplex
      #papilo
      ];
  };
  
  pyscipopt = pkgs.python3Packages.buildPythonPackage {
        pname = "PySCIPOpt";
        version = "v5.0.1";
        src = pkgs.fetchFromGitHub {
          owner = "scipopt";
          repo = "PySCIPOpt";
          rev = "v5.0.1";
          hash = "sha256-zfebn7jhvDndWAn+eB0q+2iOjPIqY143bbaKePl4FSQ=";
        };
        buildInputs = with pkgs; [ 
          python3Packages.cython
          scip 
        ];
        SCIPOPTDIR = "${scip}";
      };
      
  frobby = pkgs.stdenv.mkDerivation {
        pname = "frobby";
        version = "0.9.5";
        src = pkgs.fetchFromGitHub {
          owner = "Macaulay2";
          repo = "frobby";
          rev = "v0.9.5";
          hash = "sha256-wa1lCqQyctS3hgs84S28aRpU9UCciI0WU9oEjAzhg44=";
        };
        makeFlags = [ "PREFIX=${placeholder "out"}" "BIN_INSTALL_DIR=${placeholder "out"}/bin" ];
        postPatch = "substituteInPlace Makefile --replace '@/usr/bin/env' '@env'";
        buildInputs = with pkgs; [
          gmp
        ];
      };

in 
  pkgs.mkShell {
    buildInputs = with pkgs; [
      sage
      python3Packages.pweave
      
      z3
      python3Packages.z3
      
      frobby
      scip
      pyscipopt
    ];
}
		
