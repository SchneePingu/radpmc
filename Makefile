base=RadPMC

make:	${base}.hs
	ghc --make -O2 -dynamic -rtsopts -fforce-recomp ${base}.hs
