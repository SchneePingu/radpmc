build:  
	cabal v2-build
install:
	cabal v2-install --lib radpmc
uninstall:
	rm -rf `find ${HOME}/.cabal/ -name radpmc-*`
	sed -i '/radpmc-.*/d' `find ${HOME}/.ghc/ -name environments`/default
