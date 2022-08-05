#include <stdlib.h>
#include <time.h>
#include "stdio.h"
#include <math.h>
#include "Simulator/sim.h"
#include "Simulator/norm.h"
#include <string.h>
#include <unistd.h>
#include <stdbool.h>  

int main(int argc, char *argv[]){
    if( argc != 7 ) {
        printf("The format should be: distance(i) error(f) transveral(0/1) theta(f) phi(f) fast_mode(0/1)\n");
    }   
    
    int ancilla, nw, ne, sw, se, qubits, bit[1], logic,i, sweep, syndromeindex, xstabs, zstabs, roundlen;
    static int distance;
    sscanf(argv[1],"%d",&distance);
    double p, angleth, angleph;
    ds_Register reg;
  
    int width = 2*distance+1;

    int **grid = calloc(width, sizeof(int *));
    for (int i = 0; i < width; i++) {
        grid[i] = calloc(width, sizeof(int));
    }

    qubits = 1;
    xstabs = 0;
    zstabs = 0;
    roundlen = 0;

    for (int j = 1; j < 2*distance+1; j++){
        for (int i = 1; i < 2*distance+1; i++){
            if (i%2 && j%2){
                grid[j][i] = qubits;
                qubits++;
            }
        }
    }

    for (int j = 0; j < 2*distance+1; j++){
        for (int i = 2; i < 2*distance-1; i++){
            if ((i%4==2 && j%4==0) || (i%4==0 && j%4==2)){
                xstabs++;
                if ((i==0) || (j==0) || (i==2*distance) || (j==2*distance)) {roundlen=roundlen+4;} else {roundlen = roundlen+6;}
            }            
            if ((j%4==0 && i%4==0) || (j%4==2 && i%4==2)){
                zstabs++;
                if ((j==0) || (i==0) || (j==2*distance) || (i==2*distance)) {roundlen=roundlen+2;} else {roundlen = roundlen+4;}
            }
        }
    }

    for (int j = 2; j < 2*distance-1; j++){
        for (int i = 0; i < 2*distance+1; i++){
        }
    }

    int syndromeX[distance+1][8], syndromeZ[distance+1][7];
    sscanf(argv[2],"%le",&p);
    ds_initialize_simulator((unsigned) time(NULL) * getpid());
    reg = ds_create_register(qubits, p, 0);
    ds_set_state(reg, 0, 1, 0);



    printf("ds_create_register(%d, %f, 0);\n", qubits, p);
    printf("ds_set_state(reg, 0, 1, 0);\n");

 // Un comment for transversal injection on the bloch spehere equator (angle = rotation away from |+> state on Bloch spehere in radians)
    if (*argv[3] == '1'){
       sscanf(argv[4],"%le",&angleth);
       sscanf(argv[5],"%le",&angleph);
       for (i=1; i<=distance*distance; i++){

           ds_yrot(reg,i,angleth,0);
           ds_zrot(reg,i,angleph,0);
       }        
    }
    int runs = 0;
    bool error = false;
    double rngbuffer[roundlen*(distance+1)];

    // Fast_mode code. If enabled, will roll until an error occurs. Assumes no errors will result in a correct run. Not appropriate for post-select data.
    while (!error) {
    	for(int i = 0; i < roundlen*(distance+1); i++) {
    		if (i <= roundlen*distance) {
    			rngbuffer[i] = ds_uniform();
    		} else {
    			rngbuffer[i] = 0.999;
    		}
    		if (rngbuffer[i] < p){
    			error = true;
    		}
    	}
    	if (*argv[6] == '0') {error = true;}

    	if ((!error) && (p!=0)) {
    		runs++;
    	} else {
    		error = true;
    	}
    }

    printf("Runs: %d\n", runs);   
    bit[0] = 0;
 	int ii = 0;

    // Loop over d rounds of stabilizer measurements
    for(sweep = 0; sweep < distance + 1; sweep++) {
    	// qubit 0 is a re-used ancilla
    	// qubit n (n!=0) is the corresponding data qubit

        syndromeindex = 0;
        for (int j = 0; j < 2*distance+1; j++){
            for (int i = 2; i < 2*distance-1; i++){
                if ((i%4==2 && j%4==0) || (i%4==0 && j%4==2)){
                    ds_Hadamard(reg, 0, 1);ds_lerr(reg, 0, 1, rngbuffer[ii]);ii++;
                    if ((i-1)>0 && (j+1)<2*distance){ nw = grid[j+1][i-1]; } else { nw = -1; };
                    if ((i-1)>0 && (j-1)>0){ sw = grid[j-1][i-1]; } else { sw = -1; };
                    if ((i+1)<2*distance && (j+1)<2*distance){ ne = grid[j+1][i+1]; } else { ne = -1; };
                    if ((i+1)<2*distance && (j-1)>0){ se = grid[j-1][i+1]; } else { se = -1; };

                    if (nw>=0) {ds_cnot(reg, 0, nw, 1);ds_lerr2(reg, 0, nw, 1, rngbuffer[ii]);ii++;};
                    if (sw>=0) {ds_cnot(reg, 0, sw, 1);ds_lerr2(reg, 0, sw, 1, rngbuffer[ii]);ii++;};
                    if (ne>=0) {ds_cnot(reg, 0, ne, 1);ds_lerr2(reg, 0, ne, 1, rngbuffer[ii]);ii++;};
                    if (se>=0) {ds_cnot(reg, 0, se, 1);ds_lerr2(reg, 0, se, 1, rngbuffer[ii]);ii++;};
                    ds_Hadamard(reg, 0, 1);ds_lerr(reg, 0, 1, rngbuffer[ii]);ii++;

                    syndromeX[sweep][syndromeindex] = ds_measure(reg, 1, bit);
                    if (syndromeX[sweep][syndromeindex] == 1) ds_X(reg, 0, 1);
                    printf("ax %d,%d bit %d\n", j, i, syndromeX[sweep][syndromeindex]);
                    syndromeindex++;
                }
            }
        }
        syndromeindex = 0;
        for (int j = 2; j < 2*distance-1; j++){
            for (int i = 0; i < 2*distance+1; i++){
                if ((i%4==0 && j%4==0) || (i%4==2 && j%4==2)){
                    if ((i-1)>0 && (j+1)<2*distance){ nw = grid[j+1][i-1]; } else { nw = -1; };
                    if ((i-1)>0 && (j-1)>0){ sw = grid[j-1][i-1]; } else { sw = -1; };
                    if ((i+1)<2*distance && (j+1)<2*distance){ ne = grid[j+1][i+1]; } else { ne = -1; };
                    if ((i+1)<2*distance && (j-1)>0){ se = grid[j-1][i+1]; } else { se = -1; };

                    if (nw>=0) {ds_cnot(reg, nw, 0, 1);ds_lerr2(reg, nw, 0, 1, rngbuffer[ii]);ii++;};
                    if (sw>=0) {ds_cnot(reg, sw, 0, 1);ds_lerr2(reg, sw, 0, 1, rngbuffer[ii]);ii++;};
                    if (ne>=0) {ds_cnot(reg, ne, 0, 1);ds_lerr2(reg, ne, 0, 1, rngbuffer[ii]);ii++;};
                    if (se>=0) {ds_cnot(reg, se, 0, 1);ds_lerr2(reg, se, 0, 1, rngbuffer[ii]);ii++;};
                    syndromeZ[sweep][syndromeindex] = ds_measure(reg, 1, bit);
                    if (syndromeZ[sweep][syndromeindex] == 1) ds_X(reg, 0, 1);
                    printf("az %d,%d bit %d\n", j, i, syndromeZ[sweep][syndromeindex]);
                    syndromeindex++;
                }
            }
        }
 


        printf("Sweep: %d\n", sweep);
        // abort if error detected
        if (sweep > 0) {        	
	        for (int i = 0; i < xstabs; i++) {
	        	printf("X %d, %d, %d, %d\n", i, sweep, syndromeX[sweep-1][i], syndromeX[sweep][i]);
		        if (syndromeX[sweep][i] != syndromeX[sweep-1][i]) {
		        	exit(1);
		        }
		    }
	        for (int i = 0; i < zstabs; i++) {
	        	printf("Z %d, %d, %d, %d\n", i, sweep, syndromeZ[sweep-1][i], syndromeZ[sweep][i]);
		        if (syndromeZ[sweep][i] != syndromeZ[sweep-1][i]) {
		        	exit(1);
		        }
		    }
        }
 
    }
    
    //Print full state vector of the computer
    ds_print(reg);
    ds_destroy_register(reg);
    return 0;

}

;
/*-----------------------end-------------------------------*/

