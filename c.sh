
xcc -target=XK-EVK-XU316 -O3 -I./EmbeddedLapack/src/LinearAlgebra -I./EmbeddedLapack/src/Lapack/Include -report -c ./EmbeddedLapack/src/Lapack/Scr/dgeev.c
xcc -target=XK-EVK-XU316 -O3 -I./EmbeddedLapack/src/LinearAlgebra -I./EmbeddedLapack/src/Lapack/Include -report -c ./EmbeddedLapack/src/Lapack/Scr/dgebak.c
xcc -target=XK-EVK-XU316 -O3 -I./EmbeddedLapack/src/LinearAlgebra -I./EmbeddedLapack/src/Lapack/Include -report -c ./EmbeddedLapack/src/LinearAlgebra/eig.c 
xcc -target=XK-EVK-XU316 -std=c++14 -O3 -I./EmbeddedLapack/src/LinearAlgebra -I./EmbeddedLapack/src/Lapack/Include -report -c main.cpp
xcc -target=XK-EVK-XU316 -O3 -report -o a.xe main.o eig.o tran.o dgeev.o dgebak.o


