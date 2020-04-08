/*
Code by Gabe Shepherd
PHYS 370
Started: 16 April 2019
Completed:
Ballistics Gel Simulator
*/

// import libraries
#include <mygraph.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

// global constants
#define dim 2
#define I 60
#define J 60

// global variables
double r0[2];
double R0 = 1;
double R1;
double bpos[dim] = {5,5};
double b0pos[dim] = {-10,30};
double bvel[dim] = {0,0};
double b0vel[dim] = {1, 0};
double ppos[dim][I][J];
double pvel[dim][I][J];
double bradius = 5;
double simtime = 0;
int blacklist[8][I][J];
double Kconst = 0.01;
double forcemax = 1;
double dt = 0.001;
double pmass = 1;
double bmass = 100;
int progress = 1;
char sweepname[100] = "finalprojdata.csv";
double gscale = 20;
double CM[dim] = {0, 0};

// objects
void setCM(void){
    double MM=0;
    for( int d=0; d<dim; d++ ){
        CM[d]=0;
    }
    for( int i=0; i<I; i++ ){
        for( int j=0; j<J; j++ ){
            MM += pmass;
            for( int d=0; d<dim; d++ ){
                CM[d] += ppos[d][i][j]*pmass;
            }
        }
    }
    for( int d=0; d<dim; d++ ){
        CM[d]/=MM;
    }
}

void init(){
    R1 = sqrt(2)*R0;
    for( int i=0 ; i<8 ; i++ ){
        if( i%2 == 0 ){
            r0[i] = R0;
        }
        else if( i%2 != 0 ){
            r0[i] = R1;
        }
    }
    for( int i=0 ; i<dim ; i++ ){
        bpos[i] = b0pos[i];
        bvel[i] = b0vel[i];
    }
    for( int i=0 ; i<I ; i++ ){
        for( int j=0 ; j<J ; j++ ){
            for( int k=0 ; k<dim ; k++ ){
                if( i==0 && k==0 ){
                    ppos[k][i][j] = 0;
                }
                else if( j==0 && k==1 ){
                    ppos[k][i][j] = 0;
                }
                else{
                    if( k==0 ){
                        ppos[k][i][j] = ppos[k][i-1][j] + R0;
                    }
                    else if ( k==1 ){
                        ppos[k][i][j] = ppos[k][i][j-1] + R0;
                    }
                }
            }
        }
    }
    for( int i=0 ; i<I ; i++ ){
        for( int j=0 ; j<J ; j++ ){
            for( int k=0 ; k<dim ; k++ ){
                pvel[k][i][j] = 0;
            }
        }
    }
    for( int i=0; i<I; i++ ){
        for( int j=0; j<J; j++ ){
            for( int k=0; k<8; k++ ){
                blacklist[k][i][j] = 0;
            }
        }
    }
    for( int i=0; i<I; i++ ){
        for( int j=0; j<J; j++ ){
            for( int k=0; k<4; k++ ){
                if( i==0 ){
                    if( k==3 ){
                        blacklist[k][i][j] = 1;
                    }
                }
                if( j==J-1 ){
                    if( k==1 || k==2 || k==3 ){
                        blacklist[k][i][j] = 1;
                    }
                }
                if( i==I-1 ){
                    if( k==0 || k==1 ){
                        blacklist[k][i][j] = 1;
                    }
                }
            }
        }
    }
    setCM();
}

int iterate(){
    double momentum;
    double dist[I][J];
    int interact[I][J];
    double deltax[I][J];
    double deltay[I][J];
    double changex, changey;
    double distance;
    double deltar[8][I][J];
    double force;
    double sforce[dim][8][I][J];
    double pforce[dim][I][J];
    double bforce[dim];
    double Frepulse[I][J];
    double base;
    double deltav;
    FILE *res1=fopen("test.csv", "w");

    memset(&sforce[0][0][0][0],0,dim*8*I*J*sizeof(double));
    memset(&pforce[0][0][0],0,dim*I*J*sizeof(double));
    for (int d=0; d<dim; d++) bforce[d]=0;
    for( int i=0; i<I; i++ ){
        for( int j=0; j<J; j++ ){
            dist[i][j] = sqrt( ( (ppos[0][i][j]-bpos[0])*(ppos[0][i][j]-bpos[0]) ) + ( (ppos[1][i][j]-bpos[1])*(ppos[1][i][j]-bpos[1]) ) );
            deltax[i][j] = (ppos[0][i][j]-bpos[0]);
            deltay[i][j] = (ppos[1][i][j]-bpos[1]);
            base = dist[i][j] - bradius;
            Frepulse[i][j] = 1/pow( base, 12 );
            for( int d=0; d<dim; d++ ){
                if( d==0 ){
                    pforce[d][i][j] = Frepulse[i][j] * (deltax[i][j]/dist[i][j]);
                    bforce[d] -= pforce[d][i][j];
                }
                else if( d==1 ){
                    pforce[d][i][j] = Frepulse[i][j] * (deltay[i][j]/dist[i][j]);
                    bforce[d] -= pforce[d][i][j];
                }
            }
            for( int k=0; k<4; k++ ){
	      /* This needs to be done only once not every iteration
                if( i==0 ){
                    if( k==3 ){
                        blacklist[k][i][j] = 1;
                    }
                }
                if( j==J-1 ){
                    if( k==1 || k==2 || k==3 ){
                        blacklist[k][i][j] = 1;
                    }
                }
                if( i==I-1 ){
                    if( k==0 || k==1 ){
                        blacklist[k][i][j] = 1;
                    }
                }
		/* end comments */
                if( blacklist[k][i][j]==0 ){
                    if( k==0 ){
                        changex = ppos[0][i+1][j] - ppos[0][i][j];
                        changey = ppos[1][i+1][j] - ppos[1][i][j];
                        distance = sqrt( ( changex*changex ) + ( changey*changey ) );
                        deltar[k][i][j] = distance - r0[k];
                        force = Kconst*deltar[k][i][j];
                        for( int d=0; d<dim;d++ ){
                            if( force>forcemax ){
                                blacklist[k][i][j] = 1;
                            }
                            else if( d==0 ){
                                sforce[d][k][i][j] = force * (changex/distance);
                                sforce[d][k+4][i+1][j] = -sforce[d][k][i][j];
                            }
                            else if( d==1 ){
                                sforce[d][k][i][j] = force * (changey/distance);
                                sforce[d][k+4][i+1][j] = -sforce[d][k][i][j];
                            }
                        }
                    }
                    else if( k==1 ){
                        changex = ppos[0][i+1][j+1] - ppos[0][i][j];
                        changey = ppos[1][i+1][j+1] - ppos[1][i][j];
                        distance = sqrt( ( changex*changex ) + ( changey*changey ) );
                        deltar[k][i][j] = distance - r0[k];
                        force = Kconst*deltar[k][i][j];
                        for( int d=0; d<dim;d++ ){
                            if( force>forcemax ){
                                blacklist[k][i][j] = 1;
                            }
                            else if( d==0 ){
                                sforce[d][k][i][j] = force * (changex/distance);
                                sforce[d][k+4][i+1][j+1] = -sforce[d][k][i][j];
                            }
                            else if( d==1 ){
                                sforce[d][k][i][j] = force * (changey/distance);
                                sforce[d][k+4][i+1][j+1] = -sforce[d][k][i][j];
                            }
                        }
                    }
                    else if( k==2 ){
                        changex = ppos[0][i][j+1] - ppos[0][i][j];
                        changey = ppos[1][i][j+1] - ppos[1][i][j];
                        distance = sqrt( ( changex*changex ) + ( changey*changey ) );
                        deltar[k][i][j] = distance - r0[k];
                        force = Kconst*deltar[k][i][j];
                        for( int d=0; d<dim;d++ ){
                            if( force>forcemax ){
                                blacklist[k][i][j] = 1;
                            }
                            else if( d==0 ){
                                sforce[d][k][i][j] = force * (changex/distance);
                                sforce[d][k+4][i][j+1] = -sforce[d][k][i][j];
                            }
                            else if( d==1 ){
                                sforce[d][k][i][j] = force * (changey/distance);
                                sforce[d][k+4][i][j+1] = -sforce[d][k][i][j];
                            }
                        }
                    }
                    else if( k==3 ){
                        changex = ppos[0][i-1][j+1] - ppos[0][i][j];
                        changey = ppos[1][i-1][j+1] - ppos[1][i][j];
                        distance = sqrt( ( changex*changex ) + ( changey*changey ) );
                        deltar[k][i][j] = distance - r0[k];
                        force = Kconst*deltar[k][i][j];
                        for( int d=0; d<dim;d++ ){
                            if( force>forcemax ){
                                blacklist[k][i][j] = 1;
                            }
                            else if( d==0 ){
                                sforce[d][k][i][j] = force * (changex/distance);
                                sforce[d][k+4][i-1][j+1] = -sforce[d][k][i][j];
                            }
                            else if( d==1 ){
                                sforce[d][k][i][j] = force * (changey/distance);
                                sforce[d][k+4][i-1][j+1] = -sforce[d][k][i][j];
                            }
                        }
                    }
                }
            }
        }
    }
    for( int i=0; i<I; i++ ){
      for( int j=0; j<J; j++ ){
        for( int k=0; k<8; k++ ){
	  for( int d=0; d<dim; d++ ){
	    pforce[d][i][j] += sforce[d][k][i][j];
	  }
	  
	}
      }
    }
    for( int i=0; i<I; i++ ){
        for( int j=0; j<J; j++ ){
            for( int d=0; d<dim; d++ ){
               momentum = pforce[d][i][j];
               deltav = momentum/pmass;
               pvel[d][i][j] += deltav;
               ppos[d][i][j] += (pvel[d][i][j] * dt);
            }
        }
    }
    for( int d=0; d<dim; d++ ){
        momentum = bforce[d] * dt;
        deltav = momentum/bmass;
        bvel[d] += deltav;
        bpos[d] += (bvel[d] * dt);
    }
    setCM();
    simtime += dt;
    fclose(res1);
    return progress;
}

void draw( int xdim, int ydim ){
    for( int i=0 ; i<I ; i++ ){
        for( int j=0 ; j<J ; j++ ){
            myfilledcircle( 1, xdim/2 + (ppos[0][i][j] - CM[0])*gscale, ydim/2 - (ppos[1][i][j] - CM[1])*gscale, 10 );
        }
    }
    myfilledcircle( 2, xdim/2 + (bpos[0] - CM[0])*gscale, ydim/2 - (bpos[1] - CM[1])*gscale, 90 );
}

void main(void){
    struct timespec mytime;
    mytime.tv_sec=0;
    mytime.tv_nsec=100000;

    int done=0,cont=0,sstep=0, repeat=1;  
    AddFreedraw("Simulate", draw);
    StartMenu("Top Menu", 1);
    DefineDouble("dt", &dt);
    DefineDouble("Scale Factor", &gscale);
    DefineString("Filename",&sweepname[0],100);
    DefineGraph(freedraw_, "Simulation");
    StartMenu("Gel Info", 1);
    DefineDouble("Dist btw Particles", &R0);
    DefineDouble("Mass of Particle", &pmass);
    DefineDouble("Spring Constant", &Kconst);
    DefineDouble("Max Force", &forcemax);
    EndMenu();
    StartMenu("Bullet Info", 1);
    DefineDouble("Bullet Radius", &bradius);
    DefineDouble("Bullet Mass", &bmass);
    DefineDouble("Initial X Bullet Velocity", &b0vel[0]);
    DefineDouble("Initial Y Bullet Velocity", &b0vel[1]);
    DefineDouble("Initial X Bullet Position", &b0pos[0]);
    DefineDouble("Initial Y Bullet Position", &b0pos[1]);
    EndMenu();
    DefineFunction("Init", &init);
    DefineLong("sleep", &mytime.tv_nsec);
    DefineInt("repeat", &repeat);
    DefineBool("sstep", &sstep);
    DefineBool("cont", &cont);
    DefineBool("done", &done);
    EndMenu();

    init();
    FILE *res=fopen(sweepname, "w");
    while( !done ){
        Events(1);
        DrawGraphs();
        if( cont || sstep ){
            sstep=0;
            fprintf(res, "%e %e\n", bpos[0], bvel[0] );
            for( int i=0 ; i<repeat ; i++ ){
                progress = iterate();
                if ( progress==0 ){
                    cont = 0;
                    done = 0;
                    break;
                }
            }
            nanosleep(&mytime, NULL);
        }
    }
    fclose(res);
}
