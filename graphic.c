#include <allegro.h>
#include <math.h>
#include <stdio.h>
#include "ptask.h"
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#define W_win 640
#define H_win 480
#define X1 120
#define X2 520
#define Y1 300
#define BKG 0
#define CR_CMP_R 160
#define CR_CMP_G 82
#define CR_CMP_B 45
#define CR_All_R 105
#define CR_All_G 105
#define CR_All_B 105
#define scale 2

//valori di OFFSET bitmap MQ e P_data
#define W_mq 16
#define H_mq 7
#define X_OFF 37
#define Y_OFF 15
#define X_TEXT 35
#define Y_TEXT 5
#define X_TEXT_data 15
#define Y_TEXT_data 50
#define FB 10

//elementi per aggiornamento della tabella degli stati
#define N_ST_SV 5
#define X_Mat_S_OFF 355
#define Y_Mat_S_OFF  25
#define Y_Lab_S_OFF 5
#define X_Lab_S_OFF 405
#define L_S_rect 25
#define C_S_rect 13
#define C_S_RAD 10
#define G_X_OFF 15
#define G_Y_OFF 170
#define Len_Ax_X 300
#define Len_Ax_Y 140
//elementi da togliere perchè nella define di qlearning
#define N_state 36
#define N_state_x_ang 6
#define N_action 4

#define h_floor 50
#define w_centre 150
#define d_body 0.2
#define h_body 0.03
#define r_wheel 0.015
#define r_joint 0.0075

extern float* get_state();

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



static int figure[12];
static int Stat_lp[N_ST_SV];
static int Reward_lp[Len_Ax_X];
static int Sl_count=0;
static int Rew_count=0;
static BITMAP * CR;
static BITMAP * MQ;
static BITMAP * P_data;
static BITMAP* GRP_STAT;
static float MQ_[N_state][N_action];
int conv_col(int col,int cscale)
{
    if(col*cscale>=255)
    {
        return(255);
    }
    else
    {
        return(col*cscale);
    }
}
void   init_s()
{
    char s[20];
    int col,VB,VR,VG,c_scale=1,i,j;
    //inizializzazione finestra allegro 
    allegro_init();
    set_color_depth(32);
    set_gfx_mode(GFX_AUTODETECT_WINDOWED,W_win*scale,H_win*scale,0,0);
    clear_to_color(screen,BKG);
    CR=create_bitmap((X2-X1)*scale,Y1*scale);
    VB=255;
    VR=255;//valore da settare
    VG=conv_col(VR,c_scale);
    col=makecol(VR,VG,VB);
    clear_to_color(CR,makecol(255,255,255));
    line(CR,0,(CR->h-h_floor*scale),(CR->w),(CR->h-h_floor*scale),1);
    blit(CR,screen,0,0,X1*scale,0,CR->w,CR->h);
    MQ=create_bitmap(X1*scale,Y1*scale);
    clear_to_color(MQ,makecol(255,255,255));
    blit(MQ,screen,0,0,0,0,MQ->w,MQ->h);
    P_data=create_bitmap((W_win-X2)*scale,Y1*scale);
    clear_to_color(P_data,makecol(255,255,255));
    blit(P_data,screen,0,0,X2*scale,0,P_data->w,P_data->h);
    
    GRP_STAT=create_bitmap((W_win-X1)*scale,(Y1)*scale);
    printf("%d\n",W_win);
    clear_to_color(GRP_STAT,makecol(255,255,255));
    blit(GRP_STAT,screen,0,0,X1*scale,Y1*scale,GRP_STAT->w,GRP_STAT->h);
    //inizializzazione vettore ciclico S_loop
    for(i=0;i<N_ST_SV;i++)
    {
        Stat_lp[i]=-1;
    }
    // parte di inizializzazione pre-extern variable
    state.q1=0;
    state.q2=0;
    state.q3=0;
    state.q4=0;
    state.q5=3.14/2;
    state.q6=0;
    for(i=0;i<N_state;i++)
    {
        for(j=0;j<N_action;j++)
        {

            MQ_[i][j]=-4*i+j;
            //printf("%d ",4*i+j);
        }
        //printf("\n");
    }
    


}

//funzione che converte i metri in pixel tale che un metroo sono 1000 pixel con scala
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
void update_data(float alpha,
                float gam,
                float eps,
                float decay,
                float dis,
                int epoch ){
    char str[25];
/*
    static float alpha;		// learning rate    
    static float gam;		// discount factor
    static float epsilon; 	// actual exploration probability
    static float decay;		// decay rate for epsilon
    */
    sprintf(str,">Learning Rate:%.4f",alpha);
    textout_ex(P_data, font, str, X_TEXT_data*scale, Y_TEXT_data*scale, makecol(
0,0,0), makecol(255,255,255));
    sprintf(str,">Discount Factor:%.4f",gam);
    textout_ex(P_data, font, str, X_TEXT_data*scale, Y_TEXT_data*2*scale, makecol(
0,0,0), makecol(255,255,255));
    
    textout_ex(P_data, font,">Actual Exploration", X_TEXT_data*scale,3* Y_TEXT_data*scale, makecol(
0,0,0), makecol(255,255,255));
    sprintf(str,"Probability:%.4f",eps);
    textout_ex(P_data, font, str, X_TEXT_data*scale,(3* Y_TEXT_data+FB)*scale, makecol(
0,0,0), makecol(255,255,255));
   
    textout_ex(P_data, font,">Decay Rate for", X_TEXT_data*scale,4* Y_TEXT_data*scale, makecol(
0,0,0), makecol(255,255,255));
    sprintf(str," Epsilon:%.4f",decay);
    textout_ex(P_data, font, str,X_TEXT_data*scale,(4*Y_TEXT_data+FB)*scale, makecol(
0,0,0), makecol(255,255,255));
sprintf(str,">Distance:%.4f",dis);
    textout_ex(P_data, font, str, X_TEXT_data*scale, 5*Y_TEXT_data*scale, makecol(
0,0,0), makecol(255,255,255));
    sprintf(str,">Epoch:%d",epoch);
    textout_ex(P_data, font, str, X_TEXT_data*scale, 6*Y_TEXT_data*scale, makecol(
0,0,0), makecol(255,255,255));


blit(P_data,screen,0,0,X2*scale,0,P_data->w,P_data->h);
}


void update_STAT(int new_state)
{
    
    int i,j,k,col;
    Stat_lp[Sl_count]=new_state;
    Sl_count=(Sl_count+1)%N_ST_SV;
    textout_ex(GRP_STAT, font, "State Matrix", X_Lab_S_OFF*scale, Y_Lab_S_OFF*scale, makecol(
0,0,0), makecol(255,255,255));
    for(i=0;i<N_state_x_ang;i++)
    {
        for(j=0;j<N_state_x_ang;j++)
        {
            rect(GRP_STAT,(X_Mat_S_OFF+i*L_S_rect)*scale,
            (Y_Mat_S_OFF+j*L_S_rect)*scale,
            (X_Mat_S_OFF+(i+1)*L_S_rect)*scale,
            (Y_Mat_S_OFF+(j+1)*L_S_rect)*scale,
            makecol(0,0,0));
            for(k=0;k<N_ST_SV;k++)
            {
                //printf("Lo stato %d vale %d mentre lo stato attuale vale %d\n",k,Stat_lp[k],i*N_state_x_ang+j);
                if(Stat_lp[k]==(i*N_state_x_ang+j))
                {
                    col=245*k/N_ST_SV;
                    circlefill(GRP_STAT,
                    (X_Mat_S_OFF+i*L_S_rect+C_S_rect)*scale,
                    (Y_Mat_S_OFF+j*L_S_rect+C_S_rect)*scale,
                    C_S_RAD,
                    makecol(col,col,col)
                    );
                }
            }
            //printf("%d ",4*i+j);
        }
        //printf("\n");
    }
    
}
void update_graph(float reward,int flag,int min_range,int max_range){
   line(GRP_STAT,G_X_OFF*scale,G_Y_OFF*scale,G_X_OFF*scale,(G_Y_OFF-Len_Ax_Y)*scale,makecol(0,0,0));
   line(GRP_STAT,G_X_OFF*scale,G_Y_OFF*scale,(G_X_OFF+Len_Ax_X)*scale,G_Y_OFF*scale,makecol(0,0,0));
   if(flag==1)
   {
    if(Rew_count>Len_Ax_X)
    {
        for(int i=0;i<Len_Ax_X;i++)
        {
            Reward_lp[i]=-1;
        } 
    }
    else
    {
        if(reward>max_range)
            Reward_lp[Rew_count]=max_range;
        else if(reward<min_range)
            Reward_lp[Rew_count]=min_range;
        else
            Reward_lp[Rew_count]=reward;
        for(int i=0;i<Len_Ax_X;i++)
        {
            if(Reward_lp[i]>0)
            {
                
            }
        } 
    }
   }
}
void update_GRP_STAT(int state)
{
    clear_to_color(GRP_STAT,makecol(255,255,255));
    update_STAT(state);
    update_graph(10,1);
    blit(GRP_STAT,screen,0,0,X1*scale,Y1*scale,GRP_STAT->w,GRP_STAT->h);
    
}
//funzione di utilita per il calcolo di una posizione su un Corpo rigido 
void calc_pos(double x,double y, gsl_matrix* R,gsl_vector* p,gsl_vector*ris)
{
    gsl_vector_set(p,0,x);
    gsl_vector_set(p,1,y);
    gsl_vector_set(p,2,0);
    gsl_vector_set(p,3,1);
    gsl_blas_dgemv(CblasNoTrans,1.0,R,p,0.0,ris);
}
//funzione per il calcolo dei PX per disegnare il corpo, a seguire c'è
//quella per il link 1 e poi per il secondo
void  body_kin(int position[])
{   
    gsl_matrix* R_op;
    gsl_vector* P_op;
    gsl_vector* P_ris;
    double value;

    R_op=gsl_matrix_alloc(4,4);
    P_op=gsl_vector_alloc(4);
    P_ris=gsl_vector_alloc(4);
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
    //recupero la posizione del primo giunto direttamente dal vettore prec. calcolato
    
    position[0]=position[4];
    position[1]=position[5];

    position[2]=  MToPx(cos(state.q3)*(0.1 + 0.06 *cos(state.q4)) + 
    (-0.1 + 0.1* cos(state.q3) - 0.03*sin(state.q3)) - 
   sin(state.q3)* (0.015 + 0.06 *sin(state.q4)),0);
    position[3]= MToPx( (0.015 + 0.03 *cos(state.q3) + 0.1* sin(state.q3)) + (0.1 + 
      0.06* cos(state.q4))* sin(state.q3) + cos(state.q3) *(0.015 + 0.06* sin(state.q4)),1);
      
   /* printf("L1x:%f\n",cos(state.q3)*(0.1 + 0.06 *cos(state.q4)) + 
    (-0.1 + 0.1* cos(state.q3) - 0.03*sin(state.q3)) - 
   sin(state.q3)* (0.015 + 0.06 *sin(state.q4)));
      printf("L1y:%f\n", (0.015 + 0.03 *cos(state.q3) + 0.1* sin(state.q3)) + (0.1 + 
      0.06* cos(state.q4))* sin(state.q3) + cos(state.q3) *(0.015 + 0.06* sin(state.q4)));*/
}
void L2_kin(int position[])
{
    //recupero la posizione del primo giunto direttamente dal vettore prec. calcolato
    
    position[0]=position[2];
    position[1]=position[3];

    position[2]=  MToPx(
      (-0.1  + 0.1* cos(state.q3) - 0.03* sin(state.q3)) + 
   cos(state.q3) *(0.1 + 0.06* cos(state.q5) *sin(state.q4) + 
      cos(state.q4) *(0.06 + 0.06 *sin(state.q5))) - 
   sin(state.q3)* (0.015 - 0.06 *cos(state.q4)* cos(state.q5) + 
      sin(state.q4) *(0.06 + 0.06 *sin(state.q5)))
        ,0);
    position[3]= MToPx( 
     (0.015 + 0.03* cos(state.q3) + 0.1* sin(state.q3)) + 
   sin(state.q3)* (0.1 + 0.06* cos(state.q5)* sin(state.q4) + 
      cos(state.q4)* (0.06 + 0.06 *sin(state.q5))) + 
   cos(state.q3)* (0.015 - 0.06 *cos(state.q4)* cos(state.q5) + 
      sin(state.q4)* (0.06 + 0.06 *sin(state.q5)))
      ,1);
     /* printf("L2x:%f\n",(-0.1  + 0.1* cos(state.q3) - 0.03* sin(state.q3)) + 
   cos(state.q3) *(0.1 + 0.06* cos(state.q5) *sin(state.q4) + 
      cos(state.q4) *(0.06 + 0.06 *sin(state.q5))) - 
   sin(state.q3)* (0.015 - 0.06 *cos(state.q4)* cos(state.q5) + 
      sin(state.q4) *(0.06 + 0.06 *sin(state.q5))));
      printf("L2y:%f\n",(0.015 + 0.03* cos(state.q3) + 0.1* sin(state.q3)) + 
   sin(state.q3)* (0.1 + 0.06* cos(state.q5)* sin(state.q4) + 
      cos(state.q4)* (0.06 + 0.06 *sin(state.q5))) + 
   cos(state.q3)* (0.015 - 0.06 *cos(state.q4)* cos(state.q5) + 
      sin(state.q4)* (0.06 + 0.06 *sin(state.q5))));*/
}
void update_CR()
{   
    body_kin(figure);
    polygon(CR,5,figure,makecol(CR_CMP_R,CR_CMP_G,CR_CMP_B));
    circlefill(CR,figure[10],figure[11],MToPx(r_wheel,2),makecol(CR_All_R,CR_All_G,CR_All_B));
    circlefill(CR,figure[4],figure[5],MToPx(r_joint,2),makecol(CR_All_R,CR_All_G,CR_All_B));
    L1_kin(figure);
    line(CR,figure[0],figure[1],figure[2],figure[3],makecol(CR_CMP_R,CR_CMP_G,CR_CMP_B));
    circlefill(CR,figure[2],figure[3],MToPx(r_joint,2),makecol(CR_All_R,CR_All_G,CR_All_B));
    L2_kin(figure);
    line(CR,figure[0],figure[1],figure[2],figure[3],makecol(CR_CMP_R,CR_CMP_G,CR_CMP_B));
    blit(CR,screen,0,0,X1*scale,0,CR->w,CR->h);

}
void update_MQ(float * matrix,float step)
{
    int i,j,val;
    
    textout_ex(MQ, font, "Matrice Q", X_TEXT*scale, Y_TEXT*scale, makecol(0,0,0), makecol(255,255,255));
    for(i=0;i<N_state;i++)
    {
        for(j=0;j<N_action;j++)
        {
            if(matrix[i*N_action+j]>0)
            {
                val=(int)floor(matrix[i*N_action+j]/step);
                if(val>255)
                    val=255;
                printf("la cella %d,%d ha %d\n",i,j,val);
                rectfill(MQ,
                scale*((W_mq)*j+X_OFF),
                scale*(H_mq*i+Y_OFF),
                scale*(W_mq*(j+1)+X_OFF),
                scale*(H_mq*(i+1)+Y_OFF),
                makecol(255-val,255-val,255));
            }
            if(matrix[i*N_action+j]<0)
            {
                val=(int)floor(-matrix[i*N_action+j]/step);
                if(val>255)
                    val=255;
                printf("la cella %d,%d ha %d\n",i,j,val);
                rectfill(MQ,
                scale*(W_mq)*j+X_OFF,
                scale*H_mq*i+Y_OFF,
                scale*W_mq*(j+1)+X_OFF,
                scale*H_mq*(i+1)+Y_OFF,
                makecol(255,255-val,255-val));
            }
        }
    }
    blit(MQ,screen,0,0,0,0,MQ->w,MQ->h);
}
void update_graphic()
{
    
    
    
}
// main esistente solo per il testing 
int main()
{
    int a;
    init_s();
    update_CR();
    update_MQ((float*)&MQ_,0.5);
    update_data(0.5,0,20.2,0,0,0);
    update_GRP_STAT(10);
    update_GRP_STAT(11);
    update_GRP_STAT(12);
    update_GRP_STAT(13);
    update_GRP_STAT(14);
    update_GRP_STAT(15);
    update_GRP_STAT(16);
    update_GRP_STAT(17);
    update_GRP_STAT(18);
    update_GRP_STAT(19);
    do{
        scanf("%d",&a);
    }while(a==0); 
    printf("sono passato \n");
    allegro_exit();
    return(0);


}


