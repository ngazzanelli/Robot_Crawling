#include <allegro.h>
#include <math.h>
#include <stdio.h>
#include "ptask.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#define W_win 640
#define H_win 480
#define X1 140
#define X2 540
#define Y1 300
#define BKG 0
#define CR_CL 2
#define scale 2

#define h_floor 50
#define w_centre 150
#define d_body 0.2
#define h_body 0.03
#define r_wheel 0.015

//inserisco altre define che identificano la forma dei parametri da cercare 
struct DOF{
    float q1;
    float q2;
    float q3;
    float q4;
    float q5;
    float q6;
};

struct DOF state;

gsl_matrix* R_op;
gsl_vector* P_op;
gsl_vector* P_ris;

int figure[12];
BITMAP * CR;

void   init_s()
{
    char s[20];
    //inizializzazione finestra allegro 
    allegro_init();
    set_gfx_mode(GFX_AUTODETECT_WINDOWED,W_win*scale,H_win*scale,0,0);
    clear_to_color(screen,BKG);
    CR=create_bitmap((X2-X1)*scale,Y1*scale);
    clear_to_color(CR,15);
    line(CR,0,(CR->h-h_floor*scale),(CR->w),(CR->h-h_floor*scale),1);
    blit(CR,screen,0,0,X1*scale,0,CR->w,CR->h);
    state.q1=0;
    state.q2=0;
    state.q3=0;
    state.q4=0;
    state.q5=0;
    state.q6=0;
    R_op=gsl_matrix_alloc(4,4);
    P_op=gsl_vector_alloc(4);
    //P_int=gsl_matrix_alloc(3,1);
    P_ris=gsl_vector_alloc(4);


}
int MToPx(double val,int xy)
{
    if (xy==0)
    {
        return(((int)floor(val*1000)+w_centre)*scale);
    }
    if(xy==1)
    {
        return(CR->h-((int)floor(val*1000)+h_floor)*scale);
    }
    else{
        
        return((int)floor(val*1000*scale));
    }
}
void calc_pos(double x,double y, gsl_matrix* R,gsl_vector* p,gsl_vector*ris)
{
    gsl_vector_set(P_op,0,x);
    gsl_vector_set(P_op,1,y);
    gsl_vector_set(P_op,2,0);
    gsl_vector_set(P_op,3,1);
    gsl_blas_dgemv(CblasNoTrans,1.0,R,p,0.0,ris);
}
void  body_kin(int position[])
{   
    //cos(a),-sin(a),0
    //sin(a),cos(a),0
    //0,0,1
    double value;
    //struct fig body;
    value=cos(state.q3);
    gsl_matrix_set(R_op,0,0,value);
    value=-sin(state.q3);
    gsl_matrix_set(R_op,0,1,value);
    gsl_matrix_set(R_op,0,2,0);
    value=-0.1+0.1*cos(state.q3)-0.03*sin(state.q3);
    gsl_matrix_set(R_op,0,3,value);
    value=sin(state.q3);
    gsl_matrix_set(R_op,1,0,value);
    value=cos(state.q3);
    gsl_matrix_set(R_op,1,1,value);
    gsl_matrix_set(R_op,1,2,0);
    value=0.015+0.03*cos(state.q3)+0.1*sin(state.q3);
    gsl_matrix_set(R_op,1,3,value);
    gsl_matrix_set(R_op,2,0,0);
    gsl_matrix_set(R_op,2,1,0);
    gsl_matrix_set(R_op,2,2,1);
    gsl_matrix_set(R_op,2,3,0);
    gsl_matrix_set(R_op,3,0,0);
    gsl_matrix_set(R_op,3,1,0);
    gsl_matrix_set(R_op,3,2,0);
    gsl_matrix_set(R_op,3,3,1);

    calc_pos(-d_body/2,-h_body/2,R_op,P_op,P_ris);
    
    position[0]=MToPx(gsl_vector_get(P_ris,0),0);
    position[1]=MToPx(gsl_vector_get(P_ris,1),1);

    calc_pos(-d_body/2,h_body/2,R_op,P_op,P_ris);

    position[2]=MToPx(gsl_vector_get(P_ris,0),0);
    position[3]=MToPx(gsl_vector_get(P_ris,1),1);

    calc_pos(d_body/2,h_body/2,R_op,P_op,P_ris);

    position[4]=MToPx(gsl_vector_get(P_ris,0),0);
    position[5]=MToPx(gsl_vector_get(P_ris,1),1);

    calc_pos(d_body/2,-h_body/2-2*r_wheel,R_op,P_op,P_ris);

    position[6]=MToPx(gsl_vector_get(P_ris,0),0);
    position[7]=MToPx(gsl_vector_get(P_ris,1),1);

    calc_pos(d_body/2-r_wheel,-h_body/2,R_op,P_op,P_ris);

    position[8]=MToPx(gsl_vector_get(P_ris,0),0);
    position[9]=MToPx(gsl_vector_get(P_ris,1),1);

    calc_pos(-d_body/2,-h_body/2-r_wheel,R_op,P_op,P_ris);

    position[10]=MToPx(gsl_vector_get(P_ris,0),0);
    position[11]=MToPx(gsl_vector_get(P_ris,1),1);
    
    }
void L1_kin(int position[])
{
    
}
void update_CR()
{   
    body_kin(figure);
    polygon(CR,5,figure,CR_CL);
    circle(CR,figure[10],figure[11],MToPx(r_wheel,2),CR_CL);
    blit(CR,screen,0,0,X1*scale,0,CR->w,CR->h);

}
void update_graphic()
{
    
    
    
}
// main esistente solo per il testing 
int main()
{
    init_s();
    update_CR();
    int a;
    do{
        scanf("%d",&a);
    }while(a==0); 
    printf("sono passato \n");
    allegro_exit();
    return(0);


}

