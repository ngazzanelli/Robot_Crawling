#include <allegro.h>
#include <math.h>
#include <stdio.h>
#include "ptask.h"
#include <string.h>
#include "matrices.h"

#define W_win 640 //larghezza della finestra
#define H_win 480 //altezza della finestra 
#define W_BMS 214 //larghezza BM dei grafici 
#define H_BMS 240 //altezza BM dei grafici 
#define X_off 10 //offset a destra e sinistra del grafico 
#define Y_off 15 //offset sopra e sotto il grafico 
#define X_TEXT_OFF 107 //offset del nome del grafico 

#define PRI     10      // tasks priority 
#define DL      20      // tasks deadline
#define PER     20      // tasks period

#define BKG 0
#define num_graph 6
#define scale 2

#define PI 3.14

#define X_point 220
#define Y_point 210

//definizioni degli elementi esterni 

extern void get_sys_state(int* ic);
extern void set_sys_state(int ic);

extern void* wave_gener(void *arg);
extern void* interface(void * arg);

//qui devono essere aggiunte le extern per lo stato vero e per 
//il thread della dinamica
extern void init_state();
extern void get_state(state* rob);
extern void* dynamics(void* arg);
extern void* update_graphic(void* arg);

//makecol(255,255,255)

void print_axes(BITMAP* Graph, int G_n)
{
    char s[20];
    line(Graph,X_off*scale,(H_BMS/2)*scale,(X_point+X_off)*scale,(H_BMS/2)*scale,makecol(0,0,0));
    line(Graph,X_off*scale,Y_off*scale,X_off*scale,(Y_point+Y_off)*scale,makecol(0,0,0));
    sprintf(s,"joint var %d",G_n+1);
    textout_ex(Graph, font, s, X_TEXT_OFF*scale, Y_off*scale, makecol(0,0,0), makecol(255,255,255));
}

void init_G_DC(BITMAP* Graph[], int Graph_points[num_graph][X_point*scale], int Graph_first_p[], int Graph_elem_num[], float Graph_sym_range[])
{
    int i;
    allegro_init();
    install_keyboard();
    set_color_depth(32);
    set_gfx_mode(GFX_AUTODETECT_WINDOWED,W_win*scale,H_win*scale,0,0);
    clear_to_color(screen,BKG);

    for(i=0; i<num_graph; i++)
    {
        Graph[i]=create_bitmap(W_BMS*scale,H_BMS*scale);
        clear_to_color(Graph[i],makecol(255,255,255));
        print_axes(Graph[i],i);
        blit(Graph[i],screen,0,0,(i%3)*W_BMS*scale,(int)(i/3)*H_BMS*scale,Graph[i]->w,Graph[i]->h);
        Graph_first_p[i]=0;
        Graph_elem_num[i]=0;
        //se si volgiono mettere differenti range simmetrici
        //basta inserire degli IF/SWITCH sul numero di grafico
        Graph_sym_range[i]=PI;
    }


}

void update_gr(BITMAP* Graph,int G_n,int points[num_graph][X_point*scale],int* point_beg,int* point_num,float sym_range,float new)
{
    float  value;
    int i;
    value=new+sym_range;
    if(value>(2*sym_range))
        value=2*sym_range;
    else if(value<0)
        value=0;
    value=floor((value/(2*sym_range))*(Y_point*scale));
   // printf("il valore ridimensionato è %f\nmentre il numero originale è %f\n",value,new);

    if(*point_num<X_point*scale-1)
    {
        points[G_n][*point_num]=(int)value;
        *point_num=*point_num+1;
    }
   else
    {

        *point_beg=(*point_beg+1)%(X_point*scale);
        points[G_n][(*point_beg+X_point*scale-1)%(X_point*scale)]=(int)value;
       // printf("il valore di partenza è %d mentre l'elemento finale %d\n",Graph_first_p[G_n],(Graph_first_p[G_n]+X_point*scale)%(X_point*scale));
        //printf("il valore che è stato inserito è %d\n",Graph_points[G_n][(Graph_first_p[G_n]+(X_point*scale-1))%(X_point*scale)]);
    }
  //printf("sono passato al ciclo for dell'update\n");
    clear_to_color(Graph,makecol(255,255,255));
    print_axes(Graph,G_n);
    for(i=0;i<*point_num;i++)
    {
        putpixel(Graph,i+X_off*scale,(H_BMS-Y_off)*scale-points[G_n][(*point_beg+i)%(X_point*scale)],makecol(255,0,0));
        //printf("metto un pixel in %d, %d\n",i+X_off*scale,(H_BMS-Y_off)*scale-Graph_points[G_n][Graph_first_p[G_n]+i]);
    }
     blit(Graph,screen,0,0,(G_n%3)*W_BMS*scale,(int)(G_n/3)*H_BMS*scale,Graph->w,Graph->h);
    
    

    

}
void *update_graphic_DC(void *arg)
{
    
    int ti,exec=1,ex_stat,j;
    //float FState[6];
    BITMAP * Graph[num_graph];
    int Graph_points[num_graph][X_point*scale];
    int Graph_first_p[num_graph];
    int Graph_elem_num[num_graph];
    float Graph_sym_range[num_graph];
    printf("start graphic stack\n");
    ti = pt_get_index(arg);
    pt_set_activation(ti);
    init_G_DC(Graph,Graph_points,Graph_first_p,Graph_elem_num,Graph_sym_range);

    state rob;
    float tmp[num_graph];

    while (exec)
    {
        get_sys_state(&ex_stat);
        //if(ex_stat!=0)
         //printf("la variabile di ex_stat vale %d\n",ex_stat);
        if(ex_stat==1)
        {   //da sostituire la get per lo stato reale
            //get_FALSE_ST(FState);
            get_state(&rob);
            tmp[0] = rob.q1;
            tmp[1] = rob.q2;
            tmp[2] = rob.q3;
            tmp[3] = rob.q4;
            tmp[4] = rob.q5;
            tmp[5] = rob.q6;

            //printf("il valore di q è: [%f %f %f %f %f %f]\n",tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5]);
            for(j=0;j<num_graph;j++)
                update_gr(Graph[j],j,Graph_points,&Graph_first_p[j],&Graph_elem_num[j],Graph_sym_range[j],/*FState[j]*/tmp[j]);
            //levare il for ed utilizzare le variabili di stato vero nell'ultimo elemento della 
            //update
        }
        else if(ex_stat==3)
            exec=0;
        
        /*sezione critica per leggere le variabili di controllo 
        proposte dall'interprete e dichiarate al suo interno

        se la variabile di abilitazione del flusso è
        TRUE allora si esegue il codice altrimenti si passa direttamente
        alla condizione sulla DL_miss
        
        sezione critica per accedere alle variabili di dato,
        //in questo caso lo stato e basta.
        
        chiamata funzione per l'update della grafica*/

        pt_deadline_miss(ti);
        pt_wait_for_period(ti);
    }
    printf("end graphic task \n");
 
    return NULL;
}
/*
int main()
{
    int i,ris;
    float i_FS[6]={0,0,0,0,0,0};
    printf("inizializzo le variabili globali\n");
    // inizializzo le variabili di stato false e le var di comunicazione
    set_FALSE_ST(i_FS);
    set_com_variable(0);

    //init_state();
    printf("creo il task di interfaccia::\n");
    ris=pt_task_create( interface, 1, PER, DL, PRI);
    printf("con il risultato %d\n",ris);
    printf("creo il task di generazione d'onda\n");
    pt_task_create( wave_gener, 2, PER, DL, PRI);
    printf("creo il task di gestione della grafica\n");
    pt_task_create( update_graphic, 3, PER, DL, PRI);
    printf("creo il task per la risoluzione della dinamica\n");
    pt_task_create( dynamics, 4, 1, DL, PRI);
    for(i = 1; i <= 3; i++){
		  pt_wait_for_end(i);
		  printf("fine ciclo %d\n", i);
    }


    /*int a=0,G_n=0;
    float i;
    BITMAP * Graph[num_graph];
    int Graph_points[num_graph][X_point*scale];
    int Graph_first_p[num_graph];
    int Graph_elem_num[num_graph];
    float Graph_sym_range[num_graph];
    init_G_DC(Graph,Graph_points,Graph_first_p,Graph_elem_num,Graph_sym_range);
   for(i=-PI;i<PI;i=i+0.0025)
    {
        
        update_gr(Graph[0],0,Graph_points,&Graph_first_p[0],&Graph_elem_num[0],Graph_sym_range[0],i);
        update_gr(Graph[1],1,Graph_points,&Graph_first_p[1],&Graph_elem_num[1],Graph_sym_range[1],cos(i));
        //update_gr(,sin(i));
    }
    
    //printf("PI vale %f mentre Graph_sym_[0] vale %f\n",PI,Graph_sym_range[0]);
    do{
        scanf("%d",&a);
    }while(a==0); 
    printf("sono passato \n");
    allegro_exit();
    return 0;
}



}*/