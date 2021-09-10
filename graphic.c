#include <allegro.h>
#include <math.h>
#include <stdio.h>
#include "ptask.h"
#include <string.h>
#include "matrices.h"

// COSTANTI PER LO STATO DEL SISTEMA
#define RESET   0
#define PLAY    1
#define PAUSE   2
#define STOP    3

#define W_WIN 640 //larghezza della finestra
#define H_WIN 480 //altezza della finestra 
#define X1 120 //posizioni di riferimento delle BITMAP che compongono la parte
#define X2 520 // grafica
#define Y1 300
#define BKG 0
#define CR_CMP_R 160//colori della bitmap contenente il robot
#define CR_CMP_G 82
#define CR_CMP_B 45
#define CR_All_R 105
#define CR_All_G 105
#define CR_All_B 105
#define SCALE 2//fattore di scala

//valori di OFFSET bitmap MQ 
#define W_MQ 16 //dimensioni di ogni cella della matrice Q
#define H_MQ 5
#define X_OFF 37 //offset delle celle della matrice q
#define Y_OFF 50

#define X_TEXT 35//valori di offset della stringa "matrice Q"
#define Y_TEXT 10

#define X_TEXT_data 15//valori di offset per il plot dei 
#define Y_TEXT_data 50//dati di interesse(parametri di learning,...)
#define FB 10//offset di "a capo" qunado stampo i dati

//elementi per aggiornamento della tabella degli stati
#define N_ST_SV 5
#define X_Mat_S_OFF 370
#define Y_Mat_S_OFF  25
#define Y_Lab_S_OFF 5
#define X_Lab_S_OFF 405
#define L_S_rect 20
#define C_S_rect 10
#define C_S_RAD 10
//elementi per l'aggiornamento del grafico 
#define G_X_OFF 50
#define G_Y_OFF 165
#define Len_Ax_X 300
#define Len_Line 10
#define Len_Ax_Y 141
#define X_MAX_R_L 20
#define Y_MAX_R_L 30
#define X_MIN_R_L 20
#define Y_MIN_R_L 160
#define X_EPOCH_L 300
#define Y_EPOCH_L 170
#define X_G_NAME 100
#define Y_G_NAME 10

//elementi da togliere perchè nella define di qlearning
#define N_state 49
#define N_state_x_ang 7
#define N_action 4

//parametri per il disegno del crawler più le dimensioni fisiche
#define h_floor 50
#define w_centre 150
#define d_body 15.0
#define h_body 3.0
#define r_wheel 1.5
#define r_joint 0.75

//struttura per il passaggio della reward
typedef struct {
    int state;
    int reward;
    int flag;
} rs_for_plot;

static int graphic_dl;


//Funzioni dall'interprete
extern int get_sys_state(int* s);
extern int get_pause_graphic();
extern int get_parameter_selected();
extern void get_parameter_values(float *buff);

//Funzione dal qlearning
extern float ql_get_epsilon();

//Funzioni dal crawler
extern int angles2state(float t1, float t2);
extern void get_rs_for_plot(rs_for_plot* t);
extern void ql_get_Q(float * dest);

//Funzioni dal model
extern void get_state(state* s);


//Funzioni per gestire l'accesso al contatore delle 
//deadline di altri task
extern void get_interface_dl(int * dl_miss);
extern void get_crawler_dl(int * dl_miss);
extern void get_model_dl(int * dl_miss);

state joint_var;

int conv_col(int col,int cscale)
{
    if(col*cscale>=255)
        return(255);
    else
        return(col*cscale);
}

void reset_command()
{
    //clear_to_color(screen,makecol(200,200,200));
    rectfill(screen,0,Y1*SCALE,X1*SCALE,H_WIN*SCALE,makecol(200,200,200));
    textout_ex(screen,font,"Pulsanti di controllo:",X_TEXT_data*SCALE,(Y1+Y_TEXT +FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    textout_ex(screen,font,"E <--> Chiusura",X_TEXT_data*SCALE,(Y1+Y_TEXT+2*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    textout_ex(screen,font,"del Programma",X_TEXT_data*SCALE,(Y1+Y_TEXT+3*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    textout_ex(screen,font,"S <--> Avvio",X_TEXT_data*SCALE,(Y1+Y_TEXT+4*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    textout_ex(screen,font,"del Learning",X_TEXT_data*SCALE,(Y1+Y_TEXT+5*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    textout_ex(screen,font,"UP/DOWN <--> Cambio",X_TEXT_data*SCALE,(Y1+Y_TEXT+6*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    textout_ex(screen,font,"Par. di Apprendimento",X_TEXT_data*SCALE,(Y1+Y_TEXT+7*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    
    textout_ex(screen,font,"RIGHT <--> Incremento",X_TEXT_data*SCALE,(Y1+Y_TEXT+8*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    textout_ex(screen,font,"Par. di Apprendimento ",X_TEXT_data*SCALE,(Y1+Y_TEXT+9*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    
    textout_ex(screen,font,"LEFT <--> Decremento",X_TEXT_data*SCALE,(Y1+Y_TEXT+10*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    textout_ex(screen,font,"Par. di Apprendimento",X_TEXT_data*SCALE,(Y1+Y_TEXT+11*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));
    
}

void not_reset_command()
{
int txt_color;  
int bkg_color;  

    txt_color = makecol(0,0,255);       //blue
    bkg_color = makecol(200,200,200);   //grey

    rectfill(screen,0,Y1*SCALE,X1*SCALE,H_WIN*SCALE,makecol(200,200,200));
    textout_ex(screen,font,"Pulsanti di controllo:",X_TEXT_data*SCALE,(Y1+Y_TEXT +FB)*SCALE,txt_color,bkg_color);
    textout_ex(screen,font,"E <--> Chiusura",X_TEXT_data*SCALE,(Y1+Y_TEXT+2*FB)*SCALE, txt_color,makecol(200,200,200));
    textout_ex(screen,font,"del Programma",X_TEXT_data*SCALE,(Y1+Y_TEXT+3*FB)*SCALE, txt_color,makecol(200,200,200));
    textout_ex(screen,font,"R <--> Reset",X_TEXT_data*SCALE,(Y1+Y_TEXT+4*FB)*SCALE, txt_color,makecol(200,200,200));
    textout_ex(screen,font,"del Programma",X_TEXT_data*SCALE,(Y1+Y_TEXT+5*FB)*SCALE, txt_color,makecol(200,200,200));
    textout_ex(screen,font,"B <--> Boost",X_TEXT_data*SCALE,(Y1+Y_TEXT+6*FB)*SCALE, txt_color,makecol(200,200,200));
    //textout_ex(screen,font,"Accelleratore Learning",X_TEXT_data*SCALE,(Y1+Y_TEXT+7*FB)*SCALE, makecol(0,0,255),makecol(200,200,200)); 
    textout_ex(screen,font,"P <--> Pause/Play",X_TEXT_data*SCALE,(Y1+Y_TEXT+8*FB)*SCALE, makecol(0,0,255),makecol(200,200,200));

}
void   init_screen()
{
    //inizializzazione finestra allegro 
    allegro_init();
    install_keyboard();
    set_color_depth(32);
    set_gfx_mode(GFX_AUTODETECT_WINDOWED,W_WIN*SCALE,H_WIN*SCALE,0,0);
    clear_to_color(screen,makecol(200,200,200));
}

//funzione che converte i metri in pixel tale che un metroo sono 1000 pixel con scala

void update_data(BITMAP* BM_TXT,
                float alpha,
                float gam,
                float eps,
                float decay,
                float dis,
                int dl_mod,
                int dl_cra,
                int dl_int,
                int epoch )
{
    char str[25];
    clear_to_color(BM_TXT,makecol(255,255,255));
    sprintf(str,">Learning Rate:%.4f",alpha);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE, Y_TEXT_data*SCALE, makecol(0,0,0), makecol(255,255,255));
    sprintf(str,">Discount Factor:%.4f",gam);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE, (Y_TEXT_data+2*FB)*SCALE, makecol(0,0,0), makecol(255,255,255));
    textout_ex(BM_TXT, font,">Actual Exploration", X_TEXT_data*SCALE, (Y_TEXT_data+3*FB)*SCALE, makecol(0,0,0), makecol(255,255,255));
    sprintf(str,"Probability:%.4f",eps);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE, (Y_TEXT_data+4*FB)*SCALE, makecol(0,0,0), makecol(255,255,255));
    textout_ex(BM_TXT, font,">Decay Rate for", X_TEXT_data*SCALE, (Y_TEXT_data+5*FB)*SCALE, makecol(0,0,0), makecol(255,255,255));
    sprintf(str," Epsilon:%.4f",decay);
    textout_ex(BM_TXT, font, str,X_TEXT_data*SCALE, (Y_TEXT_data+6*FB)*SCALE, makecol(0,0,0), makecol(255,255,255));
    sprintf(str,">Distance:%.4f",dis);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE, (Y_TEXT_data+7*FB)*SCALE, makecol(0,0,0), makecol(255,255,255));
    sprintf(str,">Epoch:%d",epoch);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE,  (Y_TEXT_data+8*FB)*SCALE, makecol(0,0,0), makecol(255,255,255));

    sprintf(str,">Deadline Crawler:%d",dl_cra);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE,  (Y_TEXT_data+9*FB)*SCALE, makecol(
    0,0,0), makecol(255,255,255));
    sprintf(str,">Deadline Model:%d",dl_mod);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE,  (Y_TEXT_data+10*FB)*SCALE, makecol(
    0,0,0), makecol(255,255,255));
    sprintf(str,">Deadline Interpreter:%d",dl_int);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE,  (Y_TEXT_data+11*FB)*SCALE, makecol(
    0,0,0), makecol(255,255,255));
    sprintf(str,">Deadline Graphic:%d",graphic_dl);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE,  (Y_TEXT_data+12*FB)*SCALE, makecol(
    0,0,0), makecol(255,255,255));

    blit(BM_TXT,screen,0,0,X2*SCALE,0,BM_TXT->w,BM_TXT->h);
}

void update_data_reset(BITMAP* BM_TXT,
                float alpha,
                float gam,
                float eps_in,
                float eps_fi,
                float decay )
{
    char str[25];
    int select;
    clear_to_color(BM_TXT,makecol(255,255,255));
    select=get_parameter_selected();
    sprintf(str,">Learning Rate:%.4f",alpha);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE, Y_TEXT_data*SCALE, makecol(255,255,255), select==0? makecol(255,0,0): makecol(0,0,0));
    sprintf(str,">Discount Factor:%.4f",gam);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE, (Y_TEXT_data+1*FB)*SCALE, makecol(255,255,255), select==1? makecol(255,0,0): makecol(0,0,0));
    textout_ex(BM_TXT, font,">Decay Rate for", X_TEXT_data*SCALE, (Y_TEXT_data+2*FB)*SCALE, makecol(255,255,255), select==2? makecol(255,0,0): makecol(0,0,0));
    sprintf(str,"Epsilon:%.4f",decay);
    textout_ex(BM_TXT, font, str,X_TEXT_data*SCALE, (Y_TEXT_data+3*FB)*SCALE, makecol(255,255,255), select==2? makecol(255,0,0): makecol(0,0,0));
    sprintf(str,">Maximum Epsilon:%.4f",eps_in);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE, (Y_TEXT_data+4*FB)*SCALE, makecol(255,255,255), select==3? makecol(255,0,0): makecol(0,0,0));

    sprintf(str,">Minimum Epsilon:%.4f",eps_fi);
    textout_ex(BM_TXT, font, str, X_TEXT_data*SCALE, (Y_TEXT_data+5*FB)*SCALE, makecol(255,255,255), select==4? makecol(255,0,0): makecol(0,0,0));

    blit(BM_TXT,screen,0,0,X2*SCALE,0,BM_TXT->w,BM_TXT->h);
}

void update_STAT(BITMAP* BM_SG,int new_state,int reset)
{
    static int Sl_count=0;
    static int Sl_begin=0;
    static int Stat_lp[N_ST_SV];
    int i,j,k,col,ind;
    if(reset)
    {
        //printf("resetto le variabili di stato\n");
        Sl_count=0;
        Sl_begin=0;
    }
    //printf("GRAPHIC: il nuovo stato vale %d\n",new_state);
    //condizione in cui il vettore ha ancora almeno un elemento 
    //disponibile
    if(Sl_count<(N_ST_SV))
    {
        Stat_lp[Sl_count]=new_state;
        Sl_count++;
    }
    else//condizione in cui cominciano a ciclare gli elementi nel vettore
    {
        Stat_lp[(Sl_begin+N_ST_SV)%(N_ST_SV)]=new_state;
        Sl_begin=(Sl_begin+1)%N_ST_SV;
    }
    
    
    for(k=0;k<Sl_count;k++)
    {
        ind=(k+Sl_begin)%N_ST_SV;
        //printf("GRAPHIC: Lo stato %d vale %d\n",k,Stat_lp[ind]); 
    }
    textout_ex(BM_SG, font, "State Matrix", X_Lab_S_OFF*SCALE, Y_Lab_S_OFF*SCALE, makecol(0,0,0), makecol(255,255,255));
    for(i=0;i<N_state_x_ang;i++)
    {
        for(j=0;j<N_state_x_ang;j++)
        {
            rect(BM_SG,(X_Mat_S_OFF+i*L_S_rect)*SCALE,
            (Y_Mat_S_OFF+j*L_S_rect)*SCALE,
            (X_Mat_S_OFF+(i+1)*L_S_rect)*SCALE,
            (Y_Mat_S_OFF+(j+1)*L_S_rect)*SCALE,
            makecol(0,0,0));
            for(k=0;k<Sl_count;k++)
            {
                ind=(k+Sl_begin)%N_ST_SV;
                if(Stat_lp[ind]==(i*N_state_x_ang+j))
                {
                    col=245*(4-k)/N_ST_SV;
                    circlefill(BM_SG,
                    (X_Mat_S_OFF+i*L_S_rect+C_S_rect)*SCALE,
                    (Y_Mat_S_OFF+j*L_S_rect+C_S_rect)*SCALE,
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

void update_graph(BITMAP* BM_GS,float reward,int min_range,int max_range,int reset)
{
    static float  Reward_p[(Len_Ax_X/Len_Line)];
    static int rew_count=0;
    static int rew_begin=0;
    int i,cont_plot;
    float val;
    char s[10];
    if(reset)
    {
        //printf("resetto le variabili di stato\n");
        rew_begin=0;
        rew_count=0;
    }

    line(BM_GS,G_X_OFF*SCALE,G_Y_OFF*SCALE,G_X_OFF*SCALE,(G_Y_OFF-Len_Ax_Y)*SCALE,makecol(0,0,0));
    line(BM_GS,G_X_OFF*SCALE,G_Y_OFF*SCALE,(G_X_OFF+Len_Ax_X)*SCALE,G_Y_OFF*SCALE,makecol(0,0,0));
    sprintf(s,"%d",max_range);
    textout_ex(BM_GS, font, s, X_MAX_R_L*SCALE, Y_MAX_R_L*SCALE, makecol(0,0,0), makecol(255,255,255));
    sprintf(s,"%d",min_range);
    textout_ex(BM_GS, font, s, X_MIN_R_L*SCALE, Y_MIN_R_L*SCALE, makecol(0,0,0), makecol(255,255,255));

    textout_ex(BM_GS, font, "Epoch", X_EPOCH_L*SCALE, Y_EPOCH_L*SCALE, makecol(0,0,0), makecol(255,255,255));
    textout_ex(BM_GS, font, "Reward Plot", X_G_NAME*SCALE, Y_G_NAME*SCALE, makecol(0,0,0), makecol(255,255,255));
    if(rew_count<(Len_Ax_X/Len_Line))
    {
        Reward_p[rew_count]=reward;
        rew_count++;
    }
    else
    {
        rew_begin=(rew_begin+1)%((Len_Ax_X/Len_Line));
        Reward_p[(rew_begin+(Len_Ax_X/Len_Line))%((Len_Ax_X/Len_Line))]=reward;
    }
    for(i=0;i<rew_count;i++)
    {
        cont_plot=(i+rew_begin)%(Len_Ax_X/Len_Line);
        if(Reward_p[cont_plot]>=max_range)
            val=(float)max_range;
        else if(Reward_p[cont_plot]<=min_range)
            val=(float)min_range;
        else
            val=Reward_p[cont_plot];
        
        val=((val-min_range)/(max_range-min_range))*(Len_Ax_Y*SCALE-1);
        line(BM_GS,(G_X_OFF*SCALE+i*Len_Line*SCALE),(G_Y_OFF*SCALE-floor(val)),(G_X_OFF*SCALE+(i+1)*Len_Line*SCALE),(G_Y_OFF*SCALE-floor(val)),makecol(255,0,0));
    }

}

void update_GRP_STAT(BITMAP* BM_GS,int state,float reward,int max_r,int min_r,int flag,int reset)
{
    if(flag==1)
    {
        clear_to_color(BM_GS,makecol(255,255,255));
        update_STAT(BM_GS,state,reset);
        update_graph(BM_GS,reward,min_r,max_r,reset);
        blit(BM_GS,screen,0,0,X1*SCALE,Y1*SCALE,BM_GS->w,BM_GS->h);
    }
    if(reset)
    {
        clear_to_color(BM_GS,makecol(255,255,255));
        update_STAT(BM_GS,angles2state(0,0),reset);
        update_graph(BM_GS,0,min_r,max_r,reset);
        blit(BM_GS,screen,0,0,X1*SCALE,Y1*SCALE,BM_GS->w,BM_GS->h);
    }

}
/*
* set di funzioni per la gestione della bitmap in cui viene 
* disegnato il crawler.
* MToPx mappa le coordinate reali del modello nei pixel del disegno 
*
* ""_kin aggiornano il vettore position con i punti necessari ad allegro 
* per disegnare il crawler;  
* 
* update_crawler va ad utilizzare le precedenti funzioni per disegnare 
* effettivamente sulla bitmap il soggetto
*/
int MToPx(double val,int xy)
{
    if (xy==0)
        return(((int)round(val*10)+w_centre)*SCALE);
    if(xy==1)
        return((Y1*SCALE)-((int)round(val*10)+h_floor)*SCALE);
    else        
        return((int)round(val*10*SCALE));
}
void  body_kin(int position[],state s)
{   
    
    double rot_body[4];
    double pos_body[2];
    
    rot_body[0]=cos(s.q3);
    rot_body[1]=-sin(s.q3);
    pos_body[0]=-(15/2) + (15/2)* cos(s.q3) - 3* sin(s.q3);
    rot_body[2]=sin(s.q3);
    rot_body[3]=cos(s.q3);
    pos_body[1]=1.5 + 3* cos(s.q3) + (15/2) *sin(s.q3);
    position[0]=MToPx(rot_body[0]*(-d_body/2)+rot_body[1]*(-h_body/2)+pos_body[0],0);
    position[1]=MToPx(rot_body[2]*(-d_body/2)+rot_body[3]*(-h_body/2)+pos_body[1],1);
    position[2]=MToPx(rot_body[0]*(-d_body/2)+rot_body[1]*(h_body/2)+pos_body[0],0);
    position[3]=MToPx(rot_body[2]*(-d_body/2)+rot_body[3]*(h_body/2)+pos_body[1],1);
    position[4]=MToPx(rot_body[0]*(d_body/2)+rot_body[1]*(h_body/2)+pos_body[0],0);
    position[5]=MToPx(rot_body[2]*(d_body/2)+rot_body[3]*(h_body/2)+pos_body[1],1);
    position[6]=MToPx(rot_body[0]*(d_body/2)+rot_body[1]*(-h_body/2-2*r_wheel)+pos_body[0],0);
    position[7]=MToPx(rot_body[2]*(d_body/2)+rot_body[3]*(-h_body/2-2*r_wheel)+pos_body[1],1);
    position[8]=MToPx(rot_body[0]*(d_body/2-r_wheel)+rot_body[1]*(-h_body/2)+pos_body[0],0);
    position[9]=MToPx(rot_body[2]*(d_body/2-r_wheel)+rot_body[3]*(-h_body/2)+pos_body[1],1);
    position[10]=MToPx(rot_body[0]*(-d_body/2)+rot_body[1]*(-h_body/2-r_wheel)+pos_body[0],0);
    position[11]=MToPx(rot_body[2]*(-d_body/2)+rot_body[3]*(-h_body/2-r_wheel)+pos_body[1],1);
}
void L1_kin(int position[],state s)
{
    //recupero la posizione del primo giunto direttamente dal vettore prec. calcolato
    
    position[0]=position[4];
    position[1]=position[5];
    
    position[2]=  MToPx(
        
        -(15/2) + (15 *cos(s.q3))/2 + cos(s.q3) *(15/2 + 6 *cos(s.q4)) - 
 3.0* sin(s.q3) - sin(s.q3)*(3/2 + 6* sin(s.q4))
                ,0);
    position[3]= MToPx(1.5  + 3.0* cos(s.q3) + (15 *sin(s.q3))/2 + (15/2 + 6 *cos(s.q4)) *sin(s.q3) + 
 cos(s.q3)* (1.5 + 6 *sin(s.q4))
        ,1);
    printf("la posizione del secondo giunto vale (%f,%f)\n",-(15/2) + (15 *cos(s.q3))/2 + cos(s.q3) *(15/2 + 6 *cos(s.q4)) - 
 3.0* sin(s.q3) - sin(s.q3)*(3/2 + 6* sin(s.q4)),1.5  + 3.0* cos(s.q3) + (15 *sin(s.q3))/2 + (15/2 + 6 *cos(s.q4)) *sin(s.q3) + 
 cos(s.q3)* (1.5 + 6 *sin(s.q4)));

      
  
}
void L2_kin(int position[],state s)
{
    //recupero la posizione del primo giunto direttamente dal vettore prec. calcolato
    
    position[0]=position[2];
    position[1]=position[3];

    position[2]=  MToPx(
            -(15.0/2.0)  + (15*cos(s.q3))/2.0 - 3.0*sin(s.q3) -
			sin(s.q3)*(3.0/2.0 - 6.0*cos(s.q4 + s.q5) + 6.0*sin(s.q4)) +
			cos(s.q3)*(15.0/2.0 + 6.0*cos(s.q4) + 6.0*sin(s.q4 + s.q5))
        ,0);
    position[3]= MToPx( 
            1.5 + 3.0*cos(s.q3) + (15*sin(s.q3))/2 +
			cos(s.q3)*(3.0/2.0 - 6*cos(s.q4 + s.q5) + 6*sin(s.q4)) +
			sin(s.q3)*(15.0/2.0 + 6*cos(s.q4) + 6*sin(s.q4 + s.q5))
        ,1);
}

void update_CR(BITMAP* BM_CR,state joint_v)
{   
    //printf("SONO DENTRO UPDATE_CR\n");
    BITMAP * floor_cell=load_bitmap("floor.",NULL);
    if(floor_cell==NULL)
        printf("DIOPORCO\n");*/
    int figure[12];
    clear_to_color(BM_CR,makecol(255,255,255));
    line(BM_CR,0,(BM_CR->h-h_floor*SCALE),(BM_CR->w),(BM_CR->h-h_floor*SCALE),1);
    body_kin(figure,joint_v);
    polygon(BM_CR,5,figure,makecol(CR_CMP_R,CR_CMP_G,CR_CMP_B));
    circlefill(BM_CR,figure[10],figure[11],MToPx(r_wheel,2),makecol(CR_All_R,CR_All_G,CR_All_B));
    circlefill(BM_CR,figure[4],figure[5],MToPx(r_joint,2),makecol(CR_All_R,CR_All_G,CR_All_B));
    L1_kin(figure,joint_v);
    line(BM_CR,figure[0],figure[1],figure[2],figure[3],makecol(CR_CMP_R,CR_CMP_G,CR_CMP_B));
    circlefill(BM_CR,figure[2],figure[3],MToPx(r_joint,2),makecol(CR_All_R,CR_All_G,CR_All_B));
    L2_kin(figure,joint_v);
    line(BM_CR,figure[0],figure[1],figure[2],figure[3],makecol(CR_CMP_R,CR_CMP_G,CR_CMP_B));
    //blit(floor_cell,BM_CR,0,0,0,0,floor_cell->h,floor_cell->h);
    blit(BM_CR,screen,0,0,X1*SCALE,0,BM_CR->w,BM_CR->h);

}
/*
* Funzione che aggiorna la BITMAP contenente la Matrice Q, che dovrà essere 
* passsata come argomento della funzione sotto forma di puntatore;
* è possibile variare la gradazione dei colori con cui vengono 
* rappresetnati i valori della matrice graficamente mediante 
* l'argomento <step> della update
*/

void update_MQ(BITMAP* BM_MQ,float * matrix,float step)
{
    int i,j,val,col;
    clear_to_color(BM_MQ,makecol(255,255,255));
    textout_ex(BM_MQ, font, "Matrice Q", X_TEXT*SCALE, Y_TEXT*SCALE, makecol(0,0,0), makecol(255,255,255));
    for(i=0;i<N_state;i++)
    {
        for(j=0;j<N_action;j++)
        {
            if(matrix[i*N_action+j]>0)
            {
                    val=(int)floor(matrix[i*N_action+j]/step);
                    if(val>255)
                        val=255;
                    col=makecol(255-val,255-val,255);
            }
            else
            {
                 val=(int)floor(-matrix[i*N_action+j]/step);
                if(val>255)
                    val=255;
                col=makecol(255,255-val,255-val);
            }
            rectfill(BM_MQ,
                SCALE*((W_MQ)*j+X_OFF),
                SCALE*(H_MQ*i+Y_OFF),
                SCALE*(W_MQ*(j+1)+X_OFF),
                SCALE*(H_MQ*(i+1)+Y_OFF),
                col);
        }
    }
    blit(BM_MQ,screen,0,0,0,0,BM_MQ->w,BM_MQ->h);
}



void *update_graphic(void *arg)
{
   
    printf("GRAPHIC: task started\n");    
    int ti,int_dl,mod_dl,craw_dl,epoch = 0,exec;
    state rob;
    rs_for_plot rew_st;
    float Matrix_Q[49*4];
    float epsilon;
    float values[5];
    BITMAP *CR,*MQ,*P_data,*GRP_STAT;
    //inizializzo allegro e lo schermo 
    init_screen();
    //inizializzo le BITMAP

    CR=create_bitmap((X2-X1)*SCALE,Y1*SCALE);
    MQ=create_bitmap(X1*SCALE,Y1*SCALE);
    P_data=create_bitmap((W_WIN-X2)*SCALE,Y1*SCALE);
    GRP_STAT=create_bitmap((W_WIN-X1)*SCALE,(Y1)*SCALE);
    
    ti = pt_get_index(arg);
    pt_set_activation(ti);

    while (get_sys_state(&exec) != STOP)
    {

        if(exec==PLAY || exec==RESET/*&& !get_pause_graphic()*/) //decommenta se vuoi vedere la grafica che si ferma
        {
            //printf("DENTRO IF DI update_graphic\n");
            get_state(&rob);
            update_CR(CR,rob);
            
            get_interface_dl(&int_dl);
            get_model_dl(&mod_dl);
            get_crawler_dl(&craw_dl);
            epsilon=ql_get_epsilon();
            if(exec==PLAY)
            {
                update_data(P_data,values[0],values[1],epsilon,values[2],rob.q1,mod_dl,craw_dl,int_dl,epoch);
                not_reset_command();
            }
            if (exec==RESET)
            {
                get_parameter_values(values);
                update_data_reset(P_data,values[0],values[1],values[3],values[4],values[2]);
                ql_get_Q(Matrix_Q);
                update_MQ(MQ, Matrix_Q, 0.1);
                update_GRP_STAT(GRP_STAT,rew_st.state,rew_st.reward,50,-50,rew_st.flag,1);
                reset_command();
            }
            get_rs_for_plot(&rew_st);
            update_GRP_STAT(GRP_STAT,rew_st.state,rew_st.reward,50,-50,rew_st.flag,0);
            

            if(rew_st.flag==1)
            {
                ql_get_Q(Matrix_Q);
                update_MQ(MQ,Matrix_Q,0.1);
                epoch++;
            }
        }
        
        if(pt_deadline_miss(ti))
            graphic_dl++;
        pt_wait_for_period(ti);
        
    } 
    printf("GRAPHIC: task finished\n");
    return NULL;
}
