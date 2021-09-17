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

// Parameter selected constants
#define ALPHA	0
#define GAMMA	1
#define	DECAY	2
#define EPS_MAX	3
#define EPS_MIN	4

// Window Dimensions Constants
#define W_WIN   	640    // Window Width
#define H_WIN   	480    // Window Height 
#define X1      	120    
#define X2      	520    
#define Y1      	300    
#define BRD_THICK   3

// Scale Factor Constant 
#define SCALE   	2 

// Crawler plot Constants
#define BKG         0       
#define CR_CMP_R    14      //Crawler's body and link color (Electric Blue) 
#define CR_CMP_G    75 
#define CR_CMP_B    239
#define CR_All_R    105    //Crawler's joint and wheel color (Gray)
#define CR_All_G    105
#define CR_All_B    105
#define H_FLOOR		50     // Floor offset 
#define W_CENTRE	150    // 
#define D_BODY		15.0   // Body Widht 
#define H_BODY		3.0    // Body Height
#define R_WHEEL		1.5    // Wheel Radius 
#define R_JOINT		0.75   // Joint Radius 
#define ARM_THICK   2      // Link Thickness
#define N_BLOCKS    2       // # of big blocks shown in the floor 
#define BM_BLOCKS	4		// # of big blocks in floor_bitmap 
#define W_BLOCK     100     // length of little blocks (1 big block = 2x2 little block)
#define X			0		// flag used in metres2pixel conversion
#define Y			1		// flag used in metres2pixel conversion
#define W_LAND		600     // Background Width 
#define H_LAND		250     // Background Height 
#define TREE_SPACE	100		// distance between trees
#define TREE_RADIUS 35		// tree's crown radius
#define TREE_W		20		// tree's trunk width
#define TREE_H		100		// tree's trunk heigh
#define X_SUN       300     // Sun's x offset 
#define Y_SUN       250     // Sun's y offset 
#define RADIUS_SUN  25      // Sun's Radius 

// Text plot Constants
#define X_TEXT_DATA 15      // Text x offset
#define Y_TEXT_DATA 50      // Text y offset
#define FB 15               // Vertical offset between two string belong to different sentence 
#define NL 5                // Vertical offset between two string belong to the same sentence 

// QL_State plot Constants
#define N_ST_SV		5      // # of state show in QL_state plot 
#define X_MAT_S_OFF 370    // QL_State matrix x offset 
#define Y_MAT_S_OFF 25     // QL_state matrix y offset 
#define Y_LAB_S_OFF 5      // QL_state label x offset
#define X_LAB_S_OFF 405    // QL_state label y offset 
#define L_S_RECT    20     // Cell lenght 
#define C_S_RECT    10     // Cell Centre from vertex 
#define C_S_RAD     10     // State Circle Radius 


//Matrix Q plot Constant
#define W_MQ 		16     // Cell Widht 
#define H_MQ	 	5      // Cell Height 
#define X_OFF 		25     // Matrix Offset
#define Y_OFF 		50
#define X_TEXT 		35     //Text Offset 
#define Y_TEXT 		10

// Graphic plot Constants
#define G_X_OFF 	50      // grapich offset 
#define G_Y_OFF 	165
#define LEN_LINE    10      // width of a single reward plot 
#define LEN_AX_X 	300     // lenght of the axis 
#define LEN_AX_Y 	141
#define X_MAX_R_L 	20      //offset of max bound string 
#define Y_MAX_R_L 	30
#define X_MIN_R_L 	20      //offset of min bound string 
#define Y_MIN_R_L 	160
#define X_EPOCH_L 	300     //offset of epoch string 
#define Y_EPOCH_L 	170
#define X_G_NAME	100     //offset of graphic name 
#define Y_G_NAME 	10
#define MAX_REW     50      //max saturation bound for reward plotting
#define MIN_REW    -50      //min saturation bound for reward plotting

//elementi da togliere perchè nella define di qlearning
#define N_STATE 49
#define N_STATE_X_ANG 7
#define N_ACTION 4

//Reward and state struct to comunicate with Crawler.c
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

// Static variable for floor and landscape bitmap
static BITMAP*  floor_bitmap;
static BITMAP*	landscape_bitmap;


//---------------------------------------------
//  The Following Function initializes the
//  application's window
//---------------------------------------------
void init_screen()
{
    allegro_init();
    install_keyboard();
    set_color_depth(32);
    set_gfx_mode(GFX_AUTODETECT_WINDOWED,W_WIN*SCALE,H_WIN*SCALE,0,0);
    clear_to_color(screen,makecol(200,200,200));
}


//--------------------------------------------
// The following function draws a thick line
//--------------------------------------------
void thick_line(BITMAP *bmp, float x1, float y1, float x2, float y2, float thi, int color)
{
float   dx = x1 - x2;
float   dy = y1 - y2;
float   d = 2*sqrt(dx*dx + dy*dy);
float   ca = dx/d;
float   sa = dy/d;
int v[4*2];

    if (!d) return;

    v[0] = x1 - thi*sa; // left up
    v[1] = y1 + thi*ca;

    v[2] = x1 + thi*sa; // right up
    v[3] = y1 - thi*ca;

    v[4] = x2 + thi*sa; // right down
    v[5] = y2 - thi*ca;

    v[6] = x2 - thi*sa; // left down
    v[7] = y2 + thi*ca;
 
    polygon(bmp, 4, v, color);
}


//-------------------------------------------
//  The Following Function plot possible 
//  keyboard commands allowed in reset state   
//-------------------------------------------
void reset_command(BITMAP* bitmap_command)
{
int bkg_col;
int txt_col;
int border_col;

    bkg_col = makecol(200, 200, 200);	//grey
    txt_col = makecol(0, 0, 0);		    //black
    border_col = makecol(0,0,0);        //black

    clear_to_color(bitmap_command, border_col);
    rectfill(bitmap_command, BRD_THICK,bitmap_command->h - BRD_THICK, bitmap_command->w - BRD_THICK, BRD_THICK, bkg_col);
    clear_to_color(bitmap_command, border_col);
    rectfill(bitmap_command, BRD_THICK,bitmap_command->h - BRD_THICK, bitmap_command->w - BRD_THICK, BRD_THICK, bkg_col);
    textout_ex(bitmap_command, font, "Pulsanti di controllo:", X_TEXT_DATA*SCALE, (bitmap_command->h- Y_TEXT + FB)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "E <--> Chiusura", X_TEXT_DATA*SCALE, (Y_TEXT + 2*FB)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "del Programma", X_TEXT_DATA*SCALE, (Y_TEXT + 2*FB + NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "S <--> Avvio", X_TEXT_DATA*SCALE, (Y_TEXT + 3*FB + NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "del Learning", X_TEXT_DATA*SCALE, (Y_TEXT + 3*FB + 2*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "UP/DOWN <--> Cambio", X_TEXT_DATA*SCALE, (Y_TEXT + 4*FB + 2*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "Par. di Apprendimento", X_TEXT_DATA*SCALE, (Y_TEXT + 4*FB + 3*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "RIGHT <--> Incremento", X_TEXT_DATA*SCALE, (Y_TEXT + 5*FB + 3*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "Par. di Apprendimento ", X_TEXT_DATA*SCALE, (Y_TEXT + 5*FB + 4*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "LEFT <--> Decremento", X_TEXT_DATA*SCALE, (Y_TEXT + 6*FB + 4*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "Par. di Apprendimento", X_TEXT_DATA*SCALE, (Y_TEXT + 6*FB + 5*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "L <--> Load Q Matrix", X_TEXT_DATA*SCALE, (Y_TEXT + 7*FB + 5*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "from File", X_TEXT_DATA*SCALE, (Y_TEXT + 7*FB + 6*NL)*SCALE, txt_col, bkg_col);
    blit(bitmap_command, screen, 0, 0, 0, Y1*SCALE, bitmap_command->w, bitmap_command->h);
    
}


//-------------------------------------------
//  The Following Function plot possible 
//  keyboard commands allowed in all not 
//  reset state   
//-------------------------------------------
void not_reset_command(BITMAP* bitmap_command, int state)
{
int txt_col;  
int bkg_col;  
int border_col;
    txt_col = makecol(0, 0, 0);         //black
    bkg_col = makecol(200, 200, 200);   //grey
    border_col = makecol(0,0,0);        //black

    clear_to_color(bitmap_command, border_col);
    rectfill(bitmap_command, BRD_THICK,bitmap_command->h - BRD_THICK, bitmap_command->w - BRD_THICK, BRD_THICK, bkg_col);
    textout_ex(bitmap_command, font, "Pulsanti di controllo:", X_TEXT_DATA*SCALE, (Y_TEXT + FB)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "E <--> Chiusura", X_TEXT_DATA*SCALE, (Y_TEXT + 2*FB)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "del Programma", X_TEXT_DATA*SCALE, (Y_TEXT + 2*FB + NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "R <--> Reset", X_TEXT_DATA*SCALE, (Y_TEXT + 3*FB + NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "del Programma", X_TEXT_DATA*SCALE, (Y_TEXT + 3*FB + 2*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font, "B <--> Boost", X_TEXT_DATA*SCALE, (Y_TEXT + 4*FB + 2*NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_command, font,"P <--> Pause/Play", X_TEXT_DATA*SCALE, (Y_TEXT + 5*FB + 2*NL)*SCALE, txt_col, bkg_col);
    if(state == PAUSE){
        textout_ex(bitmap_command, font,"F <--> Save", X_TEXT_DATA*SCALE, (Y_TEXT + 6*FB + 2*NL)*SCALE, txt_col, bkg_col);
        textout_ex(bitmap_command, font,"Q Matrix to File", X_TEXT_DATA*SCALE, (Y_TEXT + 6*FB + 3*NL)*SCALE, txt_col, bkg_col);
    }
    blit(bitmap_command, screen, 0, 0, 0, Y1*SCALE, bitmap_command->w, bitmap_command->h);
}


//-------------------------------------------
//  The Following Function plots data of 
//  interest in all not-reset state  
//-------------------------------------------
void update_parameter(BITMAP* bitmap_parameter,
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
    int border_col;

    bkg_col = makecol(200, 200, 200);	//grey
    txt_col = makecol(0, 0, 0);		    //black
    border_col = makecol(0,0,0);        //black

    clear_to_color(bitmap_parameter, border_col);
    rectfill(bitmap_parameter, BRD_THICK,bitmap_parameter->h - BRD_THICK, bitmap_parameter->w - BRD_THICK, BRD_THICK, bkg_col);

    sprintf(str, ">Learning Rate:%.4f", alpha);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + FB)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Discount Factor:%.4f", gam);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 2*FB)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_parameter, font, ">Actual Exploration", X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 3*FB)*SCALE, txt_col, bkg_col);
    sprintf(str," Probability:%.4f",eps);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 3*FB+ NL)*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_parameter, font,">Decay Rate for", X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 4*FB + NL)*SCALE, txt_col, bkg_col);
    sprintf(str, " Epsilon:%.4f", decay);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 4*FB + 2*NL)*SCALE, txt_col, bkg_col);
    sprintf(str,">Distance:%.4f",dis);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 5*FB + 2*NL)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Epoch:%d", epoch);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 6*FB + 2*NL)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Deadline Crawler:%d", dl_cra);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 7*FB + 2*NL)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Deadline Model:%d", dl_mod);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 8*FB + 2*NL)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Deadline Interpreter:%d", dl_int);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 9*FB + 2*NL)*SCALE, txt_col, bkg_col);
    sprintf(str, ">Deadline Graphic:%d", dl_gra);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE,  (Y_TEXT_DATA + 10*FB + 2*NL)*SCALE, txt_col, bkg_col);

    blit(bitmap_parameter, screen, 0, 0, X2*SCALE, 0, bitmap_parameter->w, bitmap_parameter->h);
}


//-------------------------------------------
//  The Following Function plots data of 
//  interest in reset state and show which 
//  is going to be modified by the user 
//-------------------------------------------
void update_parameter_reset(BITMAP* bitmap_parameter,
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
    int border_col;

    bkg_col = makecol(200, 200, 200);	//grey
    txt_col = makecol(0, 0, 0);			//black
    slc_col = makecol(0, 0, 0);			//red
    border_col = makecol(0,0,0);    	//black
	//Border drawing
    clear_to_color(bitmap_parameter, border_col);
    rectfill(bitmap_parameter, BRD_THICK,bitmap_parameter->h - BRD_THICK, bitmap_parameter->w - BRD_THICK, BRD_THICK, bkg_col);
	//Data variables printing
    select = get_parameter_selected();
    sprintf(str, ">Learning Rate: %.2f", alpha);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, Y_TEXT_DATA*SCALE, select == ALPHA ? bkg_col : txt_col, select == ALPHA ? slc_col : bkg_col);
    sprintf(str, ">Discount Factor: %.2f", gam);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 1*FB)*SCALE, select == GAMMA ? bkg_col : txt_col, select == GAMMA ? slc_col : bkg_col);
    textout_ex(bitmap_parameter, font,">Decay Rate for", X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 2*FB)*SCALE, select == DECAY ? bkg_col : txt_col, select == DECAY ? slc_col : bkg_col);
    sprintf(str, "Epsilon: %.2f", decay);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 2*FB + NL)*SCALE, select == DECAY ? bkg_col : txt_col, select == DECAY ? slc_col : bkg_col);
    sprintf(str,">Maximum Epsilon: %.2f",eps_in);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 3*FB + NL)*SCALE, select == EPS_MAX ? bkg_col : txt_col, select == EPS_MAX ? slc_col : bkg_col);
    sprintf(str, ">Minimum Epsilon: %.2f", eps_fi);
    textout_ex(bitmap_parameter, font, str, X_TEXT_DATA*SCALE, (Y_TEXT_DATA + 4*FB + NL)*SCALE, select == EPS_MIN ? bkg_col : txt_col, select == EPS_MIN ? slc_col : bkg_col);

    blit(bitmap_parameter, screen, 0, 0, X2*SCALE, 0, bitmap_parameter->w, bitmap_parameter->h);
}


//-------------------------------------------
//  The Following Function updates the plot
//  of the ql_state memorizing the last 
//  N_ST_SV ql_state
//-------------------------------------------
void update_state(BITMAP* bitmap_rew_state, int new_state, int reset)
{
    static int sl_count = 0;
    static int sl_begin = 0;
    static int stat_lp[N_ST_SV];
    int i, j, k, col, ind, txt_col, bkg_col;

    bkg_col = makecol(255, 255, 255);	//grey
    txt_col = makecol(0, 0, 0);			//white

    /*-------------------------------------------
      Circular Array update adding a new state 
      if it isn't full or replacing the 
      oldest state memorized in it 
    ------------------------------------------*/
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
    
    /*-------------------------------------------
      Plot N_STATE Cell and if the 
      corrispond state is memorize into the 
      array plot a circle inside, the level
     of brightness corrensond to the oldness 
     of the state 
    -------------------------------------------*/
    for(k = 0; k < sl_count; k++)
        ind = (k + sl_begin) % N_ST_SV; 

    textout_ex(bitmap_rew_state, font, "State Matrix", X_LAB_S_OFF*SCALE, Y_LAB_S_OFF*SCALE, txt_col, bkg_col);
    for(i = 0; i < N_STATE_X_ANG; i++){
        for(j = 0; j < N_STATE_X_ANG; j++){
            rect(bitmap_rew_state, (X_MAT_S_OFF + i*L_S_RECT)*SCALE,
                        (Y_MAT_S_OFF+j*L_S_RECT)*SCALE,
                        (X_MAT_S_OFF+(i+1)*L_S_RECT)*SCALE,
                        (Y_MAT_S_OFF+(j+1)*L_S_RECT)*SCALE,
                        txt_col
			);

            for(k = 0; k < sl_count ; k++){
                ind = (k + sl_begin) % N_ST_SV;
                if(stat_lp[ind] == (i * N_STATE_X_ANG + j)){
                    col = 245 * (4 - k) / N_ST_SV;
                    circlefill(bitmap_rew_state, 
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
void update_graph(BITMAP* bitmap_rews_state, float reward, int min_range, int max_range, int reset)
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

    line(bitmap_rews_state,G_X_OFF*SCALE, G_Y_OFF*SCALE, G_X_OFF*SCALE, (G_Y_OFF-LEN_AX_Y)*SCALE, ax_col);
    line(bitmap_rews_state,G_X_OFF*SCALE, G_Y_OFF*SCALE, (G_X_OFF+LEN_AX_X)*SCALE, G_Y_OFF*SCALE, ax_col);
    sprintf(s, "%d", max_range);
    textout_ex(bitmap_rews_state, font, s, X_MAX_R_L*SCALE, Y_MAX_R_L*SCALE, txt_col,bkg_col);
    sprintf(s,"%d", min_range);
    textout_ex(bitmap_rews_state, font, s, X_MIN_R_L*SCALE, Y_MIN_R_L*SCALE, txt_col, bkg_col);

    textout_ex(bitmap_rews_state, font, "Epoch", X_EPOCH_L*SCALE, Y_EPOCH_L*SCALE, txt_col, bkg_col);
    textout_ex(bitmap_rews_state, font, "Reward Plot", X_G_NAME*SCALE, Y_G_NAME*SCALE, txt_col, bkg_col);


    /*-------------------------------------------
      Circular Array update adding a new reward 
      if it isn't full or replacing the 
      oldest reward memorized in it 
    -------------------------------------------*/
    if(rew_count < (LEN_AX_X/LEN_LINE)){
        reward_p[rew_count] = reward;
        rew_count++;
    }else{
        rew_begin = (rew_begin + 1) % ((LEN_AX_X/LEN_LINE));
        reward_p[(rew_begin + (LEN_AX_X/LEN_LINE)) % ((LEN_AX_X/LEN_LINE))] = reward;
    }

    /*--------------------------------------------
      Plot all the reward memorized in the 
      circular array rescale it according to 
      the lenght of Y axis and the saturation 
      bound 
    --------------------------------------------*/
    for(i = 0; i < rew_count; i++){
        cont_plot = (i + rew_begin) % (LEN_AX_X/LEN_LINE);
        if(reward_p[cont_plot] >= max_range)
            val = (float)max_range;
        else if(reward_p[cont_plot] <= min_range)
            val = (float)min_range;
        else
            val = reward_p[cont_plot];
        
        val = ((val - min_range)/(max_range - min_range))*(LEN_AX_Y*SCALE - 1);
        line(bitmap_rews_state,(G_X_OFF*SCALE + i*LEN_LINE*SCALE), (G_Y_OFF*SCALE - floor(val)),
			(G_X_OFF*SCALE + (i + 1)*LEN_LINE*SCALE), (G_Y_OFF*SCALE - floor(val)), plot_col);
    }
    if(min_range <= 0 && max_range >= 0){
        val =(float)(0 - min_range)/((max_range - min_range))*(LEN_AX_Y*SCALE - 1);
        line(bitmap_rews_state, (G_X_OFF*SCALE), (G_Y_OFF*SCALE - floor(val)),
            (G_X_OFF + LEN_AX_X)*SCALE, (G_Y_OFF*SCALE - floor(val)), ax_col);
    }
}


//-------------------------------------------
//  The Following Function manages the 
//  update of reward plot and joints state
//  matrix in reset and not reset mode
//-------------------------------------------
void update_rew_state(BITMAP* bitmap_rews_state, int state, float reward, int max_r, int min_r, int reset)
{
    int border_col;
    int bkg_col;

	border_col = makecol(0,0,0);    //black
	bkg_col = makecol(255,255,255); //white
	clear_to_color(bitmap_rews_state, border_col);
	rectfill(bitmap_rews_state, BRD_THICK, bitmap_rews_state->h - BRD_THICK, bitmap_rews_state->w - BRD_THICK, BRD_THICK, bkg_col);
	if(reset){
		update_state(bitmap_rews_state, angles2state(0,0), reset);
		update_graph(bitmap_rews_state, 0, min_r, max_r, reset);
	}else{
		update_state(bitmap_rews_state, state, reset);
		update_graph(bitmap_rews_state, reward, min_r, max_r, reset);
	}
    blit(bitmap_rews_state, screen, 0, 0, X1*SCALE, Y1*SCALE, bitmap_rews_state->w, bitmap_rews_state->h);
    
}


//-------------------------------------------
//  The Following Function converts the 
//  size of the crawler part from cm to 
//  pixel
//-------------------------------------------
int cm2px(double val, int flag)
{
    if (flag == X)
        return(((int)round(val*10) + W_CENTRE)*SCALE);
    if (flag == Y)
        return((Y1*SCALE) - ((int)round(val*10) + H_FLOOR)*SCALE);
    else        
        return((int)round(val*10*SCALE));
}


//-------------------------------------------
//  The Following Functions produces the 
//  point in Pixel to draw the crawler from 
//  its kinematics
//-------------------------------------------
void  body_kin(int position[], state s)
{   
	double rot_body[4];
	double pos_body[2];

    rot_body[0] = cos(s.q3);
    rot_body[1] = -sin(s.q3);
    pos_body[0] = -(15/2) + (15/2)* cos(s.q3) - 3* sin(s.q3);
    rot_body[2] = sin(s.q3);
    rot_body[3] = cos(s.q3);
    pos_body[1] = 1.5 + 3* cos(s.q3) + (15/2) *sin(s.q3);
    position[0] = cm2px(rot_body[0]*(-D_BODY/2) + rot_body[1]*(-H_BODY/2) + pos_body[0], 0);
    position[1] = cm2px(rot_body[2]*(-D_BODY/2) + rot_body[3]*(-H_BODY/2) + pos_body[1], 1);
    position[2] = cm2px(rot_body[0]*(-D_BODY/2) + rot_body[1]*(H_BODY/2) + pos_body[0], 0);
    position[3] = cm2px(rot_body[2]*(-D_BODY/2) + rot_body[3]*(H_BODY/2) + pos_body[1], 1);
    position[4] = cm2px(rot_body[0]*(D_BODY/2) + rot_body[1]*(H_BODY/2) + pos_body[0], 0);
    position[5] = cm2px(rot_body[2]*(D_BODY/2) + rot_body[3]*(H_BODY/2) + pos_body[1], 1);
    position[6] = cm2px(rot_body[0]*(D_BODY/2) + rot_body[1]*(-H_BODY/2-2*R_WHEEL) + pos_body[0], 0);
    position[7] = cm2px(rot_body[2]*(D_BODY/2) + rot_body[3]*(-H_BODY/2-2*R_WHEEL) + pos_body[1], 1);
    position[8] = cm2px(rot_body[0]*(D_BODY/2-R_WHEEL) + rot_body[1]*(-H_BODY/2) + pos_body[0], 0);
    position[9] = cm2px(rot_body[2]*(D_BODY/2-R_WHEEL) + rot_body[3]*(-H_BODY/2) + pos_body[1], 1);
    position[10] = cm2px(rot_body[0]*(-D_BODY/2) + rot_body[1]*(-H_BODY/2-R_WHEEL) + pos_body[0], 0);
    position[11] = cm2px(rot_body[2]*(-D_BODY/2) + rot_body[3]*(-H_BODY/2-R_WHEEL) + pos_body[1], 1);
}

void L1_kin(int position[],state s)
{
    position[0] = position[4];
    position[1] = position[5];
    position[2] = cm2px(-(15/2) + (15*cos(s.q3))/2 + cos(s.q3)*(15/2 + 6*cos(s.q4)) - 3.0*sin(s.q3) - sin(s.q3)*(3/2 + 6*sin(s.q4)), 0);
    position[3] = cm2px(1.5  + 3.0*cos(s.q3) + (15*sin(s.q3))/2 + (15/2 + 6*cos(s.q4))*sin(s.q3) + cos(s.q3)*(1.5 + 6*sin(s.q4)), 1);  
}

void L2_kin(int position[],state s)
{
    position[0] = position[2];
    position[1] = position[3];
    position[2] = cm2px(
            -(15.0/2.0)  + (15*cos(s.q3))/2.0 - 3.0*sin(s.q3) -
			sin(s.q3)*(3.0/2.0 - 6.0*cos(s.q4 + s.q5) + 6.0*sin(s.q4)) +
			cos(s.q3)*(15.0/2.0 + 6.0*cos(s.q4) + 6.0*sin(s.q4 + s.q5)),0);
    position[3] = cm2px( 
            1.5 + 3.0*cos(s.q3) + (15*sin(s.q3))/2 +
			cos(s.q3)*(3.0/2.0 - 6*cos(s.q4 + s.q5) + 6*sin(s.q4)) +
			sin(s.q3)*(15.0/2.0 + 6*cos(s.q4) + 6*sin(s.q4 + s.q5)),1);
}


//-------------------------------------------
//  The Following Function updates the 
//  drawing of crawler
//-------------------------------------------
void update_crawler(BITMAP* bitmap_craw, state joint_v)
{   
int figure[12];
int x_floor_offset, x_land_offset, sun_col, sun_border,border_col, body_col, wheel_col;
    sun_col = makecol(255, 255, 0);         //yellow
    sun_border = makecol(255, 165, 0);      //orange
    border_col = makecol(0,0,0);            //black
    body_col = makecol(CR_CMP_R, CR_CMP_G, CR_CMP_B);
    wheel_col = makecol(CR_All_R, CR_All_G, CR_All_B);


    // Landscape drawing
	x_land_offset = cm2px(joint_v.q1/10, 2) % (TREE_SPACE*SCALE) + TREE_SPACE*SCALE;
	blit(landscape_bitmap, bitmap_craw, x_land_offset, 0, 0, 0, bitmap_craw->w, landscape_bitmap->h);

    // Sun drawing
    circlefill(bitmap_craw, X_SUN*SCALE,(Y1 - Y_SUN)*SCALE, RADIUS_SUN*SCALE, sun_col);
    circle(bitmap_craw, X_SUN*SCALE, (Y1 - Y_SUN)*SCALE, RADIUS_SUN*SCALE, sun_border);

    // Floor drawing
   	x_floor_offset = cm2px(joint_v.q1, 2)%(W_BLOCK*2*SCALE)+W_BLOCK*2*SCALE;
	blit(floor_bitmap, bitmap_craw, x_floor_offset, 0, 0, bitmap_craw->h - H_FLOOR*SCALE, bitmap_craw->w, floor_bitmap->h);
    line(bitmap_craw, 0, (bitmap_craw->h - H_FLOOR*SCALE), (bitmap_craw->w), (bitmap_craw->h - H_FLOOR*SCALE), 1);

    // Body drawing
    body_kin(figure, joint_v);
    polygon(bitmap_craw, 5, figure, body_col);
    line(bitmap_craw, figure[0], figure[1], figure[2], figure[3], border_col);
    line(bitmap_craw, figure[2], figure[3], figure[4], figure[5], border_col);
    line(bitmap_craw, figure[4], figure[5], figure[6], figure[7], border_col);
    line(bitmap_craw, figure[6], figure[7], figure[8], figure[9], border_col);
    line(bitmap_craw, figure[8], figure[9], figure[0], figure[1], border_col);

    // Wheel drawing
    circlefill(bitmap_craw, figure[10], figure[11], cm2px(R_WHEEL, 2), border_col); 
	circlefill(bitmap_craw, figure[10], figure[11], cm2px(R_WHEEL, 2)/2, wheel_col);

    // First link drawing
    circlefill(bitmap_craw, figure[4], figure[5], cm2px(R_JOINT ,2), wheel_col);
    circle(bitmap_craw, figure[4], figure[5], cm2px(R_JOINT, 2), border_col);
    L1_kin(figure, joint_v);
    thick_line(bitmap_craw, figure[0], figure[1], figure[2], figure[3], ARM_THICK*SCALE, body_col);

    // Second link drawing
    circlefill(bitmap_craw, figure[2], figure[3], cm2px(R_JOINT, 2), wheel_col);
    circle(bitmap_craw, figure[2], figure[3], cm2px(R_JOINT, 2), border_col);
    L2_kin(figure, joint_v);
    thick_line(bitmap_craw, figure[0],figure[1],figure[2],figure[3],ARM_THICK*SCALE,body_col);

	// Border drawing
    rectfill(bitmap_craw, 0, bitmap_craw->h,BRD_THICK,0,border_col);
    rectfill(bitmap_craw, 0, bitmap_craw->h,bitmap_craw->w,bitmap_craw->h - BRD_THICK, border_col);
    rectfill(bitmap_craw, 0, BRD_THICK,bitmap_craw->w, 0, border_col);
    rectfill(bitmap_craw, bitmap_craw->w - BRD_THICK, bitmap_craw-> h, bitmap_craw->w, 0, border_col);
	
    blit(bitmap_craw, screen, 0, 0, X1*SCALE, 0, bitmap_craw->w, bitmap_craw->h );

}

//-----------------------------------------
//  The Following Function updates the plot
//  of the Q Matrix, the parameter step 
//  allow to tune the intensity of color 
//  according on the values into Q Matrix
//-----------------------------------------
void update_qmatrix(BITMAP* bitmap_qmatrix, float* matrix, float step)
{
	int i, j, val, col, txt_col, bkg_col, border_col;
    
    bkg_col = makecol(200,200,200); //light grey
    txt_col = makecol(0,0,0);       //black
    border_col = makecol(0,0,0);    //black
    
    clear_to_color(bitmap_qmatrix, border_col);
    rectfill(bitmap_qmatrix, BRD_THICK, bitmap_qmatrix->h - BRD_THICK, bitmap_qmatrix->w - BRD_THICK, BRD_THICK, bkg_col);
    textout_ex(bitmap_qmatrix, font, "Matrice Q", X_TEXT*SCALE, Y_TEXT*SCALE, txt_col, bkg_col);

    for(i = 0; i < N_STATE; i++){
        for(j = 0; j < N_ACTION; j++){

            if(matrix[i*N_ACTION + j] > 0){
                    val = (int)floor(sqrt(matrix[i*N_ACTION + j])/step);
                    if(val > 255)
                        val=255;
                    col = makecol(255 - val, 255 - val, 255);
            }
            else{
                val = (int)floor(-matrix[i*N_ACTION + j]/step);
                if(val > 255)
                    val = 255;
                col = makecol(255, 255 - val, 255 - val);
            }
            rectfill(bitmap_qmatrix,
                SCALE*((W_MQ)*j + X_OFF),
                SCALE*(H_MQ*i + Y_OFF),
                SCALE*(W_MQ*(j + 1) + X_OFF),
                SCALE*(H_MQ*(i + 1) + Y_OFF),
                col);
            rect(bitmap_qmatrix,
                SCALE*((W_MQ)*j + X_OFF),
                SCALE*(H_MQ*i + Y_OFF),
                SCALE*(W_MQ*(j + 1) + X_OFF),
                SCALE*(H_MQ*(i + 1) + Y_OFF),
                txt_col);
        }
    }
    blit(bitmap_qmatrix, screen, 0, 0, 0, 0, bitmap_qmatrix->w, bitmap_qmatrix->h);
}

//-----------------------------------------
//	The Following Function initializes the
//	floor bitmap with predefined pattern
//-----------------------------------------
void init_floor_bitmap()
{
    int i, x1, x2;
    int floor_color1, floor_color2;

    floor_color1 = makecol(100, 200, 100);	//light green
    floor_color2 = makecol(100, 250, 100);	//dark green
	clear_to_color(floor_bitmap, floor_color1);
    for(i = 0; i < BM_BLOCKS; i++){
		x1 = (i*2*W_BLOCK + W_BLOCK)*SCALE;
		x2 = i*2*W_BLOCK*SCALE;
		//printf("GRAPHIC: la x1 vale %d mentre x2 vale %d\n", x1, x2);
        rectfill(floor_bitmap, x1, H_FLOOR/2*SCALE, ((i + 1)*2*W_BLOCK)*SCALE, 0, floor_color2);
		rectfill(floor_bitmap, x2, H_FLOOR*SCALE, ((i + 1)*2*W_BLOCK - W_BLOCK)*SCALE, H_FLOOR/2*SCALE, floor_color2);
    }
}

//-----------------------------------------
//	The Following Function draws a tree 
//	in the landscape bitmap
//-----------------------------------------
void draw_tree(int i)
{
	int x1, x2, y1, y2, xc, yc;
	int trunk_col = makecol(91, 58, 41);	//dark brown
	int crown_col = makecol(49, 127, 67);	//dark green
    int edge_col = makecol(0,0,0);          //black
	x1 = (i*TREE_SPACE-TREE_W/2)*SCALE;
	x2 = (i*TREE_SPACE+TREE_W/2)*SCALE;
	y1 = H_LAND*SCALE;
	y2 = (H_LAND - TREE_H)*SCALE;
	xc = i*TREE_SPACE*SCALE;
	yc = y2 - TREE_RADIUS*SCALE;

	rectfill(landscape_bitmap, x1, y1, x2, y2, trunk_col );
    rect(landscape_bitmap, x1, y1, x2, y2, edge_col );
	circlefill(landscape_bitmap, xc, yc, TREE_RADIUS*SCALE, crown_col);
    circle(landscape_bitmap, xc, yc, TREE_RADIUS*SCALE, edge_col);
}

//-----------------------------------------
//	The Following Function initializes the
//	floor bitmap with predefined pattern
//-----------------------------------------
void init_landscape_bitmap()
{
    int i;
	int sky;

	sky = makecol(8, 232, 222); //light blue
	clear_to_color(landscape_bitmap, sky);
    for(i = 0; i < 7; i++)
			draw_tree(i);
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
    BITMAP *crawler_bitmap, *qmatrix_bitmap, *parameter_bitmap, *rew_state_bitmap, *command_bitmap;

    //inizializzo allegro e lo schermo 
    init_screen();
    //inizializzo le BITMAP
    crawler_bitmap = create_bitmap((X2 - X1)*SCALE, Y1*SCALE);
    qmatrix_bitmap = create_bitmap(X1*SCALE, Y1*SCALE);
    parameter_bitmap = create_bitmap((W_WIN - X2)*SCALE, Y1*SCALE);
    rew_state_bitmap = create_bitmap((W_WIN - X1)*SCALE, Y1*SCALE);
    command_bitmap = create_bitmap(X1*SCALE,(H_WIN-Y1)*SCALE);
    floor_bitmap = create_bitmap(2*W_BLOCK*BM_BLOCKS*SCALE, H_FLOOR*SCALE);
	landscape_bitmap = create_bitmap(W_LAND*SCALE, H_LAND*SCALE);
	//printf("GRAPHIC: le dimensioni della floor bitmap sono %d, %d\n", floor_bitmap->h, floor_bitmap->w);
    init_floor_bitmap();
	init_landscape_bitmap();

    ti = pt_get_index(arg);
    pt_set_activation(ti);

    while (get_sys_state(&exec) != STOP){

        if(exec != PAUSE) {
            get_rs_for_plot(&rew_st);

            if(rew_st.flag == 1)
                epoch++;

            if(!get_pause_graphic()){
                get_state(&rob);
                update_crawler(crawler_bitmap,rob);
                
                int_dmiss = pt_get_dmiss(INTERFACE);
                mod_dmiss = pt_get_dmiss(MODEL);
                craw_dmiss = pt_get_dmiss(CRAWLER);
                grap_dmiss = pt_get_dmiss(GRAPHIC);
                epsilon = ql_get_epsilon();

                if(exec == PLAY){
                    update_parameter(parameter_bitmap, values[0], values[1], epsilon,values[2], rob.q1, mod_dmiss, craw_dmiss, int_dmiss, grap_dmiss, epoch);
                    not_reset_command(command_bitmap, PLAY);

					if(rew_st.flag){ 	               
                		update_rew_state(rew_state_bitmap, rew_st.state, rew_st.reward, MAX_REW, MIN_REW, 0);
                    	ql_get_Q(Matrix_Q);
                    	update_qmatrix(qmatrix_bitmap,Matrix_Q, 0.1);
                	}
                    
                }else{ //(exec == RESET)
                    get_parameter_values(values);
                    update_parameter_reset(parameter_bitmap, values[0], values[1], values[3], values[4], values[2]);
                    ql_get_Q(Matrix_Q);
                    update_qmatrix(qmatrix_bitmap, Matrix_Q, 0.1);
	                update_rew_state(rew_state_bitmap, rew_st.state, rew_st.reward, MAX_REW, MIN_REW, 1);
                    reset_command(command_bitmap);
                }
            }
        } 
        else
            not_reset_command(command_bitmap, PAUSE);
        
        pt_deadline_miss(ti);
        pt_wait_for_period(ti);
        
    } 
    printf("GRAPHIC: task finished\n");
    return NULL;
}
