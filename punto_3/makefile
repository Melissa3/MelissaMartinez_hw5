parametros.pdf: graficas.py parametros.csv
	python graficas.py

parametros.csv: a.out
	./a.out

a.out: lotka_volterra.c lotka_volterra_obs.csv
	gcc lotka_volterra.c -lm -o a.out

clean:
	rm a.out 
