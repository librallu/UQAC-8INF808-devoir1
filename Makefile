
exe: PSO.cpp Entete.h
	g++ $^ -o $@

run:
	./exe 10 0.1 0.1 0.1 100

