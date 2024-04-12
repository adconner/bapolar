let
  pkgs = import (fetchTarball "https://github.com/NixOS/nixpkgs/tarball/nixos-23.11") { config = {}; overlays = []; };
in
  pkgs.mkShell {
    buildInputs = with pkgs; [
      sage
      python3Packages.pweave
      
      (stdenv.mkDerivation (final: {
        pname = "frobby";
        version = "0.9.5";
        src = fetchFromGitHub {
          owner = "Macaulay2";
          repo = "frobby";
          rev = "v${final.version}";
          hash = "sha256-wa1lCqQyctS3hgs84S28aRpU9UCciI0WU9oEjAzhg44=";
        };
        makeFlags = [ "PREFIX=${placeholder "out"}" "BIN_INSTALL_DIR=${placeholder "out"}/bin" ];
        postPatch = "substituteInPlace Makefile --replace '@/usr/bin/env' '@env'";
        buildInputs = [
          gmp
        ];
      }))
      
    ];
}
		
