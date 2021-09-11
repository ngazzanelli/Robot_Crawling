#include <allegro.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ptask.h"
#include "matrices.h"

// System State constnts
#define RESET   	0
#define PLAY    	1
#define PAUSE   	2
#define STOP    	3

// Task Constants
#define INTERFACE   1
#define GRAPHIC     2
#define CRAWLER     3
#define MODEL       4

// Window Dimensions Constants
#define W_WIN   	640 
#define H_WIN   	480 
#define X1      	120 
#define X2      	520 
#define Y1      	300

// Scale Factor Constant 
#define SCALE   	2 

// Crawler plot Constants
#define BKG         0
#define CR_CMP_R    160
#define CR_CMP_G    82
#define CR_CMP_B    45
#define CR_All_R    105
#define CR_All_G    105
#define CR_All_B    105
#define H_FLOOR		50
#define W_CENTRE	150
#define D_BODY		15.0
#define H_BODY		3.0
#define R_WHEEL		1.5
#define R_JOINT		0.75

// Text plot Constants
#define X_TEXT_DATA 15
#define Y_TEXT_DATA 50
#define FB 10

// QL_State plot Constants
#define N_ST_SV		5
#define X_MAT_S_OFF 370
#define Y_MAT_S_OFF 25
#define Y_LAB_S_OFF 5
#define X_LAB_S_OFF 405
#define W_MQ 		16 
#define H_MQ	 	5
#define X_OFF 		37 
#define Y_OFF 		50
#define X_TEXT 		35
#define Y_TEXT 		10

// Graphic plot Constants
#define G_X_OFF 	50
#define G_Y_OFF 	165
#define LEN_AX_X 	300
#define LEN_LINE 	10
#define LEN_AX_Y 	141
#define X_MAX_R_L 	20
#define Y_MAX_R_L 	30
#define X_MIN_R_L 	20
#define Y_MIN_R_L 	160
#define X_EPOCH_L 	300
#define Y_EPOCH_L 	170
#define X_G_NAME	100
#define Y_G_NAME 	10
#define L_S_RECT 	20
#define C_S_RECT	10
#define C_S_RAD 	10

//elementi da togliere perchè nella define di qlearning
#define N_STATE 49
#define N_STATE_X_ANG 7
#define N_ACTION 4

//Reward struct to comunicate with Crawler.c
typedef struct {
    int state;
    int reward;
    int flag;
} rs_for_plot;

//Extern Functions from command_interface
extern int get_sys_state(int* s);
extern int get_pause_graphic();
extern int get_parameter_selected();
extern void get_parameter_values(float *buff);

//Extern Functions from crawler
extern int angles2state(float t1, float t2);
extern void get_rs_for_plot(rs_for_plot* t);
extern void ql_get_Q(float* dest);
extern float ql_get_epsilon();

//Extern Function from model
extern void get_state(state* s);

//Local variable to save crawler state
state joint_var;

//-------------------------------------------
// 
//-------------------------------------------
int conv_col(int col, int cscale)
{
    if(col*cscale >= 255)
        return(255);
    else
        return(col*cscale);
}

//-------------------------------------------
//  The Following Function plot possible 
//  keyboard commands allowed in reset state   
//-------------------------------------------
void reset_command()
{
int bkg_col;
int txt_col;

    bkg_col = makecol(200, 200, 200);	//grey
    txt_col = makecol(0, 0, 255);		//blue
    rectfill(screen, 0, Y1*SCALE, X1*SCALE, H_WIN*SCALE, bkg_col);
    textout_ex(screen, font, "Pulsanti di controllo:", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "E <--> Chiusura", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 2*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "del Programma", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 3*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "S <--> Avvio", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 4*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "del Learning", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 5*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "UP/DOWN <--> Cambio", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 6*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "Par. di Apprendimento", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 7*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "RIGHT <--> Incremento", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 8*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "Par. di Apprendimento ", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 9*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "LEFT <--> Decremento", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 10*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "Par. di Apprendimento", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 11*FB)*SCALE, txt_col, bkg_col);
    
}


//-------------------------------------------
//  The Following Function plot possible 
//  keyboard commands allowed in all not 
//  reset state   
//-------------------------------------------
void not_reset_command()
{
int txt_col;  
int bkg_col;  

    txt_col = makecol(0, 0, 255);       //blue
    bkg_col = makecol(200, 200, 200);   //grey

    rectfill(screen, 0, Y1*SCALE, X1*SCALE, H_WIN*SCALE, bkg_col);
    textout_ex(screen, font, "Pulsanti di controllo:", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "E <--> Chiusura", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 2*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "del Programma", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 3*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "R <--> Reset", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 4*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "del Programma", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 5*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font, "B <--> Boost", X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 6*FB)*SCALE, txt_col, bkg_col);
    textout_ex(screen, font,"P <--> Pause/Play",X_TEXT_DATA*SCALE, (Y1 + Y_TEXT + 7*FB)*SCALE, txt_col, bkg_col);

}

//---------------------------------------------
//  The Following Function initializes the
//  application's window
//---------------------------------------------
void   init_screen()
{
    allegro_init();
    install_keyboard();
    set_color_depth(32);
    set_gfx_mode(GFX_AUTODETECT_WINDOWED,W_WIN*SCALE,H_WIN*SCALE,0,0);
    clear_to_color(screen,makecol(200,200,200));
}


//-------------------------------------------
//  The Following Function plots data of 
//  interest in all not-reset state  
//-------------------------------------------
void update_data(BITMAP* BM_TXT,
                float alpha,
                float gam,
                float eps,
                float decay,
                float dis,
                int dl_mod,
                int dl_cra,
                int dl_int,
                int dl_gra,
                int epoch )
{
    char str[25];
    int bkg_col;
    int txt_col;

    bkg_col = makecol(200, 200, 200);	//grey
    txt_col = makecol(0, 0, 255);		//blue

    clear_to_color(BM_TXT, bkg_col);
    sprintf(str, ">Learning Rate:%.4f", alpha);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, Y_TEXT_DATA*SCALE, txt_col, bkg_col);
    sprintf(str, ">Discount Factor:%.4f", gam);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 2*FB)*SCALE, txt_col, bkg_col);
    textout_ex(BM_TXT, font, ">Actual Exploration", X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 3*FB)*SCALE, txt_col, bkg_col);
    sprintf(str,"Probability:%.4f",eps);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 4*FB)*SCALE, txt_col, bkg_col);
    textout_ex(BM_TXT, font,">Decay Rate for", X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 5*FB)*SCALE, txt_col, bkg_col);
    sprintf(str, " Epsilon:%.4f", decay);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 6*FB)*SCALE, txt_col, bkg_col);
    sprintf(str,">Distance:%.4f",dis);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 7*FB)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Epoch:%d", epoch);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 8*FB)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Deadline Crawler:%d", dl_cra);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 9*FB)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Deadline Model:%d", dl_mod);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 10*FB)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Deadline Interpreter:%d", dl_int);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 11*FB)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Deadline Graphic:%d", dl_gra);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE,  (Y_TEXT_DATA + 12*FB)*SCALE, txt_col, bkg_col);

    blit(BM_TXT, screen, 0, 0, X2*SCALE, 0, BM_TXT->w, BM_TXT->h);
}

//-------------------------------------------
//  The Following Function plots data of 
//  interest in reset state and show which 
//  is going to be modified by the user 
//-------------------------------------------
void update_data_reset(BITMAP* BM_TXT,
                float alpha,
                float gam,
                float eps_in,
                float eps_fi,
                float decay )
{
    
    char str[25];
    int select;
    int txt_col;  
    int bkg_col;
    int slc_col;

    bkg_col = makecol(200, 200, 200);	//grey
    txt_col = makecol(0, 0, 255);		//blue
    slc_col = makecol(255, 0, 0);		//red

    clear_to_color(BM_TXT, bkg_col);
    select = get_parameter_selected();
    sprintf(str, ">Learning Rate:%.4f", alpha);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, Y_TEXT_DATA*SCALE, txt_col, select == 0 ? slc_col : bkg_col);
    sprintf(str, ">Discount Factor:%.4f", gam);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 1*FB)*SCALE, txt_col, select == 1 ? slc_col : bkg_col);
    textout_ex(BM_TXT, font,">Decay Rate for", X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 2*FB)*SCALE, txt_col, select == 2 ? slc_col : bkg_col);
    sprintf(str, "Epsilon:%.4f", decay);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 3*FB)*SCALE, txt_col, select == 2 ? slc_col : bkg_col);
    sprintf(str,">Maximum Epsilon:%.4f",eps_in);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 4*FB)*SCALE, txt_col, select == 3 ? slc_col : bkg_col);
    sprintf(str, ">Minimum Epsilon:%.4f", eps_fi);
    textout_ex(BM_TXT, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 5*FB)*SCALE, txt_col, select == 4 ? slc_col : bkg_col);

    blit(BM_TXT, screen, 0, 0, X2*SCALE, 0, BM_TXT->w, BM_TXT->h);
}


//-------------------------------------------
//  The Following Function updates the plot
//  of the ql_state memorizing the last 
//  N_ST_SV ql_state
//-------------------------------------------
void update_STAT(BITMAP* BM_SG, int new_state, int reset)
{
    static int sl_count = 0;
    static int sl_begin = 0;
    static int stat_lp[N_ST_SV];
    int i, j, k, col, ind, txt_col, bkg_col;

    bkg_col = makecol(255, 255, 255);	//grey
    txt_col = makecol(0, 0, 0);			//white

    if(reset){
        sl_count = 0;
        sl_begin = 0;
    }
    if(sl_count < (N_ST_SV))
    {
        stat_lp[sl_count] = new_state;
        sl_count++;
    }else{
        stat_lp[(sl_begin + N_ST_SV) % (N_ST_SV)] = new_state;
        sl_begin = (sl_begin + 1) % N_ST_SV;
    }
    
    
    for(k = 0; k < sl_count; k++)
        ind = (k + sl_begin) % N_ST_SV; 

    textout_ex(BM_SG, font, "State Matrix", X_LAB_S_OFF*SCALE, Y_LAB_S_OFF*SCALE, txt_col, bkg_col);
    for(i = 0; i < N_STATE_X_ANG; i++){
        for(j = 0; j < N_STATE_X_ANG; j++){
            rect(BM_SG, (X_MAT_S_OFF + i*L_S_RECT)*SCALE,
                        (Y_MAT_S_OFF+j*L_S_RECT)*SCALE,
                        (X_MAT_S_OFF+(i+1)*L_S_RECT)*SCALE,
                        (Y_MAT_S_OFF+(j+1)*L_S_RECT)*SCALE,
                        txt_col
			);

            for(k = 0; k < sl_count ; k++){
                ind = (k + sl_begin) % N_ST_SV;
                if(stat_lp[ind] == (i * N_STATE_X_ANG + j)){
                    col = 245 * (4 - k) / N_ST_SV;
                    circlefill(BM_SG, 
						(X_MAT_S_OFF + i*L_S_RECT + C_S_RECT)*SCALE,
                    	(Y_MAT_S_OFF + j*L_S_RECT + C_S_RECT)*SCALE,
                    	C_S_RAD, makecol(col,col,col)
					);
                }
            }   
        }
    }   
}


//-------------------------------------------
//  The Following Function updates the 
//  reward plot
//-------------------------------------------
void update_graph(BITMAP* BM_GS, float reward, int min_range, int max_range, int reset)
{
    static float	reward_p[(LEN_AX_X / LEN_LINE)];
    static int		rew_count =0, rew_begin = 0;
    int				i, cont_plot, txt_col, bkg_col, ax_col, plot_col;
    float			val;
    char			s[10];

    bkg_col = makecol(255, 255, 255);
    ax_col = txt_col = makecol(0, 0, 0);
    plot_col = makecol(255, 0, 0);
    if(reset){
        rew_begin = 0;
        rew_count = 0;
    }

    line(BM_GS,G_X_OFF*SCALE,G_Y_OFF*SCALE,G_X_OFF*SCALE,(G_Y_OFF-LEN_AX_Y)*SCALE,ax_col);
    line(BM_GS,G_X_OFF*SCALE,G_Y_OFF*SCALE,(G_X_OFF+LEN_AX_X)*SCALE,G_Y_OFF*SCALE,ax_col);
    sprintf(s,"%d",max_range);
    textout_ex(BM_GS, font, s, X_MAX_R_L*SCALE, Y_MAX_R_L*SCALE, txt_col,bkg_col);
    sprintf(s,"%d",min_range);
    textout_ex(BM_GS, font, s, X_MIN_R_L*SCALE, Y_MIN_R_L*SCALE, txt_col, bkg_col);

    textout_ex(BM_GS, font, "Epoch", X_EPOCH_L*SCALE, Y_EPOCH_L*SCALE, txt_col, bkg_col);
    textout_ex(BM_GS, font, "Reward Plot", X_G_NAME*SCALE, Y_G_NAME*SCALE, txt_col, bkg_col);
    if(rew_count < (LEN_AX_X/LEN_LINE)){
        reward_p[rew_count] = reward;
        rew_count++;
    }else{
        rew_begin = (rew_begin+1)%((LEN_AX_X/LEN_LINE));
        reward_p[(rew_begin+(LEN_AX_X/LEN_LINE))%((LEN_AX_X/LEN_LINE))] = reward;
    }
    for(i=0;i<rew_count;i++){
        cont_plot = (i+rew_begin)%(LEN_AX_X/LEN_LINE);
        if(reward_p[cont_plot]>=max_range)
            val = (float)max_range;
        else if(reward_p[cont_plot]<=min_range)
            val = (float)min_range;
        else
            val = reward_p[cont_plot];
        
        val = ((val-min_range)/(max_range-min_range))*(LEN_AX_Y*SCALE-1);
        line(BM_GS,(G_X_OFF*SCALE+i*LEN_LINE*SCALE), (G_Y_OFF*SCALE-floor(val)),
			(G_X_OFF*SCALE+(i+1)*LEN_LINE*SCALE),(G_Y_OFF*SCALE-floor(val)),plot_col);
    }

}


//-------------------------------------------
//  The Following Function manages the 
//  update of reward plot and ql_state 
//  plot in reset and not reset mode
//-------------------------------------------
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

//-------------------------------------------
//  The Following Function converts the 
//  size of the crawler part from cm to 
//  Pixel
//-------------------------------------------
int MToPx(double val,int xy)
{
    if (xy==0)
        return(((int)round(val*10)+W_CENTRE)*SCALE);
    if (xy==1)
        return((Y1*SCALE)-((int)round(val*10)+H_FLOOR)*SCALE);
    else        
        return((int)round(val*10*SCALE));
}

//-------------------------------------------
//  The Following Functions produces the 
//  point in Pixel to draw the crawler from 
//  its kinematics
//-------------------------------------------
void  body_kin(int position[],state s)
{   
	double rot_body[4];
	double pos_body[2];

    rot_body[0] = cos(s.q3);
    rot_body[1] = -sin(s.q3);
    pos_body[0] = -(15/2) + (15/2)* cos(s.q3) - 3* sin(s.q3);
    rot_body[2] = sin(s.q3);
    rot_body[3] = cos(s.q3);
    pos_body[1] = 1.5 + 3* cos(s.q3) + (15/2) *sin(s.q3);
    position[0] = MToPx(rot_body[0]*(-D_BODY/2)+rot_body[1]*(-H_BODY/2)+pos_body[0],0);
    position[1] = MToPx(rot_body[2]*(-D_BODY/2)+rot_body[3]*(-H_BODY/2)+pos_body[1],1);
    position[2] = MToPx(rot_body[0]*(-D_BODY/2)+rot_body[1]*(H_BODY/2)+pos_body[0],0);
    position[3] = MToPx(rot_body[2]*(-D_BODY/2)+rot_body[3]*(H_BODY/2)+pos_body[1],1);
    position[4] = MToPx(rot_body[0]*(D_BODY/2)+rot_body[1]*(H_BODY/2)+pos_body[0],0);
    position[5] = MToPx(rot_body[2]*(D_BODY/2)+rot_body[3]*(H_BODY/2)+pos_body[1],1);
    position[6] = MToPx(rot_body[0]*(D_BODY/2)+rot_body[1]*(-H_BODY/2-2*R_WHEEL)+pos_body[0],0);
    position[7] = MToPx(rot_body[2]*(D_BODY/2)+rot_body[3]*(-H_BODY/2-2*R_WHEEL)+pos_body[1],1);
    position[8] = MToPx(rot_body[0]*(D_BODY/2-R_WHEEL)+rot_body[1]*(-H_BODY/2)+pos_body[0],0);
    position[9] = MToPx(rot_body[2]*(D_BODY/2-R_WHEEL)+rot_body[3]*(-H_BODY/2)+pos_body[1],1);
    position[10] = MToPx(rot_body[0]*(-D_BODY/2)+rot_body[1]*(-H_BODY/2-R_WHEEL)+pos_body[0],0);
    position[11] = MToPx(rot_body[2]*(-D_BODY/2)+rot_body[3]*(-H_BODY/2-R_WHEEL)+pos_body[1],1);
}
void L1_kin(int position[],state s)
{
    position[0] = position[4];
    position[1] = position[5];
    position[2] = MToPx(-(15/2) + (15 *cos(s.q3))/2 + cos(s.q3) *(15/2 + 6 *cos(s.q4)) - 3.0* sin(s.q3) - sin(s.q3)*(3/2 + 6* sin(s.q4)),0);
    position[3] = MToPx(1.5  + 3.0* cos(s.q3) + (15 *sin(s.q3))/2 + (15/2 + 6 *cos(s.q4)) *sin(s.q3) + cos(s.q3)* (1.5 + 6 *sin(s.q4)),1);  
}
void L2_kin(int position[],state s)
{
    position[0] = position[2];
    position[1] = position[3];
    position[2] = MToPx(
            -(15.0/2.0)  + (15*cos(s.q3))/2.0 - 3.0*sin(s.q3) -
			sin(s.q3)*(3.0/2.0 - 6.0*cos(s.q4 + s.q5) + 6.0*sin(s.q4)) +
			cos(s.q3)*(15.0/2.0 + 6.0*cos(s.q4) + 6.0*sin(s.q4 + s.q5)),0);
    position[3] = MToPx( 
            1.5 + 3.0*cos(s.q3) + (15*sin(s.q3))/2 +
			cos(s.q3)*(3.0/2.0 - 6*cos(s.q4 + s.q5) + 6*sin(s.q4)) +
			sin(s.q3)*(15.0/2.0 + 6*cos(s.q4) + 6*sin(s.q4 + s.q5)),1);
}

//-------------------------------------------
//  The Following Function updates the 
//  drawing of crawler
//-------------------------------------------
void update_CR(BITMAP* BM_CR,state joint_v)
{   
	int figure[12];
    //printf("SONO DENTRO UPDATE_CR\n");
	/*
	BITMAP * floor_cell;
    floor_cell=load_bitmap("sky_sfodo.bmp",NULL);
    if(floor_cell==NULL)
        printf("DIOPORCO\n");
	*/
    
    clear_to_color(BM_CR,makecol(255,255,255));
    line(BM_CR,0,(BM_CR->h-H_FLOOR*SCALE),(BM_CR->w),(BM_CR->h-H_FLOOR*SCALE),1);
    body_kin(figure,joint_v);
    polygon(BM_CR,5,figure,makecol(CR_CMP_R,CR_CMP_G,CR_CMP_B));
    circlefill(BM_CR,figure[10],figure[11],MToPx(R_WHEEL,2),makecol(CR_All_R,CR_All_G,CR_All_B));
    circlefill(BM_CR,figure[4],figure[5],MToPx(R_JOINT,2),makecol(CR_All_R,CR_All_G,CR_All_B));
    L1_kin(figure,joint_v);
    line(BM_CR,figure[0],figure[1],figure[2],figure[3],makecol(CR_CMP_R,CR_CMP_G,CR_CMP_B));
    circlefill(BM_CR,figure[2],figure[3],MToPx(R_JOINT,2),makecol(CR_All_R,CR_All_G,CR_All_B));
    L2_kin(figure,joint_v);
    line(BM_CR,figure[0],figure[1],figure[2],figure[3],makecol(CR_CMP_R,CR_CMP_G,CR_CMP_B));
    //blit(floor_cell,BM_CR,0,0,0,0,floor_cell->w,floor_cell->h);
    blit(BM_CR,screen,0,0,X1*SCALE,0,BM_CR->w,BM_CR->h);

}

//-----------------------------------------
//  The Following Function updates the plot
//  of the Q Matrix, the parameter step 
//  allow to tune the intensity of color 
//  according on the values into Q Matrix
//-----------------------------------------
void update_MQ(BITMAP* BM_MQ,float * matrix,float step)
{
	int i,j,val,col,txt_col,bkg_col;
    
    bkg_col = makecol(255,255,255);
    txt_col = makecol(0,0,0);
    clear_to_color(BM_MQ,bkg_col);
    textout_ex(BM_MQ, font, "Matrice Q", X_TEXT*SCALE, Y_TEXT*SCALE, txt_col, bkg_col);
    for(i=0;i<N_STATE;i++)
    {
        for(j=0;j<N_ACTION;j++)
        {
            if(matrix[i*N_ACTION+j]>0)
            {
                    val = (int)floor(matrix[i*N_ACTION+j]/step);
                    if(val>255)
                        val=255;
                    col = makecol(255-val,255-val,255);
            }
            else
            {
                val = (int)floor(-matrix[i*N_ACTION+j]/step);
                if(val>255)
                    val = 255;
                col = makecol(255,255-val,255-val);
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


//--------------------------------------------
// Graphic Task
//-------------------------------------------
void *update_graphic(void *arg)
{
   
    printf("GRAPHIC: task started\n");    
    int ti, int_dmiss, mod_dmiss, craw_dmiss, grap_dmiss, epoch = 0, exec;
    state rob;
    rs_for_plot rew_st;
    float Matrix_Q[49*4];
    float epsilon;
    float values[5];
    BITMAP *CR, *MQ, *P_data, *GRP_STAT;

    //inizializzo allegro e lo schermo 
    init_screen();
    //inizializzo le BITMAP

    CR = create_bitmap((X2 - X1)*SCALE, Y1*SCALE);
    MQ = create_bitmap(X1*SCALE, Y1*SCALE);
    P_data = create_bitmap((W_WIN - X2)*SCALE, Y1*SCALE);
    GRP_STAT = create_bitmap((W_WIN - X1)*SCALE, Y1*SCALE);
    
    ti = pt_get_index(arg);
    pt_set_activation(ti);

    while (get_sys_state(&exec) != STOP)
    {
        if(exec==PLAY || exec==RESET/*&& !get_pause_graphic()*/) //decommenta se vuoi vedere la grafica che si ferma
        {
            //printf("GRAPHIC: dentro if di update_graphic\n");
            get_state(&rob);
            update_CR(CR,rob);
            
            int_dmiss = pt_get_dmiss(INTERFACE);
            mod_dmiss = pt_get_dmiss(MODEL);
            craw_dmiss = pt_get_dmiss(CRAWLER);
            grap_dmiss = pt_get_dmiss(GRAPHIC);
            epsilon = ql_get_epsilon();

            if(exec == PLAY)
            {
                update_data(P_data, values[0], values[1], epsilon,values[2], rob.q1, mod_dmiss, craw_dmiss, int_dmiss, grap_dmiss, epoch);
                not_reset_command();
            }

            if (exec == RESET)
            {
                get_parameter_values(values);
                update_data_reset(P_data, values[0], values[1], values[3], values[4], values[2]);
                ql_get_Q(Matrix_Q);
                update_MQ(MQ, Matrix_Q, 0.1);
                update_GRP_STAT(GRP_STAT, rew_st.state, rew_st.reward, 50, -50, rew_st.flag, 1);
                reset_command();
            }

            get_rs_for_plot(&rew_st);
            update_GRP_STAT(GRP_STAT, rew_st.state, rew_st.reward, 50, -50, rew_st.flag, 0);
            
            if(rew_st.flag == 1)
            {
                ql_get_Q(Matrix_Q);
                update_MQ(MQ,Matrix_Q, 0.1);
                epoch++;
            }
        }
        
        pt_deadline_miss(ti);
        pt_wait_for_period(ti);
        
    } 
    printf("GRAPHIC: task finished\n");
    return NULL;
}
