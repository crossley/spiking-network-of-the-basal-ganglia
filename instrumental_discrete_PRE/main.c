//
//  main.c
//  instrumental_discrete_PRE
//
//  Created by Matthew Crossley on 3/18/14.
//  Copyright (c) 2014 Matthew Crossley. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

# pragma mark -- Function Declarations --
void simulate_acquisition_full(int simulation);
void simulate_acquisition_partial(int simulation);
void simulate_extinction(int simulation);
void initialize_parameters(void);
void allocate_buffers(void);
void initialize_buffers(void);
void record_trial(int simulation, int trial);
void reset_trial(void);
void reset_simulation(void);
void reset_records(void);
void compute_average_results(void);
void save_average_essentials(void);

void update_pf(int tick);
void update_tan(int tick);
void update_sensory(int tick);
void update_msn(int tick);
void update_motor(int tick);
void update_response(int tick);
void update_reward_full(int trial);
void update_reward_partial(int trial);
void update_reward_extinction(int trial);
void update_predicted_reward(int trial);
void update_DA(int trial);
void update_w_ctx_msn(int trial);
void update_w_pf_tan(int trial);

float cap(float val);
float pos(float val);
float randn(float mean, float variance);

# pragma mark -- Parameter Declarations --
int num_simulations;
int num_trials;
int num_steps;

int condition_flag;

int trial;
int simulation;

int num_acquisition_trials;
int num_extinction_trials;

int cue_onset;
int cue_duration;

# pragma mark -- Spiking Network Declarations --
float tau;
float noise;

float *spike;
int spike_length;
float alpha_func_a;
float alpha_func_b;

float sensory_amp;
float *sensory_square;

float *msn_v;
float *msn_u;
float *msn_output;
int *msn_spikes;

float *motor_v;
float *motor_output;
int *motor_spikes;

float pause_mod_amp;
float pause_decay;
float* pause_mod;

float pf_cue_duration;
float pf_amp;
float *pf_square;

float *tan_v;
float *tan_u;
float *tan_output;
int *tan_spikes;

float w_pf_tan_init;
float w_ctx_msn_init;
float w_pf_tan;
float w_tan_msn;
float w_ctx_msn;
float w_msn_mot;

# pragma mark -- Response and Learning Declarations --
int response;
int response_time;
int reward;
float predicted_reward;
float pr_alpha;
float rpe;
float DA;

float response_threshold;
float AMPA_threshold;
float NMDA_threshold;
float LTP_msn;
float LTD_msn;
float LTP_tan;
float LTD_tan;
float DA_base;
float w_delta;

float alpha_func_a_camkII;
float alpha_func_b_camkII;
float spike_length_camkII;
float *spike_camkII;
float *msn_camkII;
float *tan_camkII;

// record buffers - each trial of each simulation is held separately
int **response_record;
int **response_time_record;
int **reward_record;
float **predicted_reward_record;
float **rpe_record;
float **DA_record;
float **w_ctx_msn_record;
float **w_pf_tan_record;

float ***sensory_square_record;
float ***msn_v_record;
float ***msn_output_record;
int ***msn_spikes_record;
float ***pf_square_record;
float ***tan_v_record;
float ***tan_output_record;
int ***tan_spikes_record;
float ***motor_v_record;
float ***motor_output_record;
int ***motor_spikes_record;

// record buffers - average over simulations
float *response_record_ave;
float *response_time_record_ave;
float *reward_record_ave;
float *predicted_reward_record_ave;
float *rpe_record_ave;
float *DA_record_ave;
float *w_ctx_msn_record_ave;
float *w_pf_tan_record_ave;

float **sensory_square_record_ave;
float **msn_v_record_ave;
float **msn_output_record_ave;
float **msn_spikes_record_ave;
float **pf_square_record_ave;
float **tan_v_record_ave;
float **tan_output_record_ave;
float **tan_spikes_record_ave;
float **motor_v_record_ave;
float **motor_output_record_ave;
float **motor_spikes_record_ave;

# pragma mark -- MAIN --
int main(int argc, const char * argv[]) {

    printf("Beginning Simulation\n");
    
    initialize_parameters();
    allocate_buffers();
    initialize_buffers();
    
    // simulate continuous reinforcement condition
    reset_records();
    condition_flag = 0; simulation = 0; trial = 0;
    for(int j=0; j<num_simulations; j++) {
        reset_simulation();
        simulate_acquisition_full(j);
        simulate_extinction(j);
        simulation++;
    }
    compute_average_results();
    save_average_essentials();
    
    // simulate partial reinforcement condition
    reset_records();
    condition_flag = 1; simulation = 0; trial = 0;
    for(int j=0; j<num_simulations; j++) {
        reset_simulation();
        simulate_acquisition_partial(j);
        simulate_extinction(j);
        simulation++;
    }
    compute_average_results();
    save_average_essentials();
    
    printf("Finished Simulation\n");
    
    return 0;
}

void simulate_acquisition_full(int simulation) {
    int jj = simulation;
    for(int i=0; i<num_acquisition_trials; i++) {
        if(trial % 10 == 0) printf("Simulating trial %i\n", trial);
        reset_trial();
        for(int j=1; j<num_steps; j++) {
            update_pf(j);
            update_tan(j);
            update_sensory(j);
            update_msn(j);
            update_motor(j);
            update_response(j);
//            if(response) break;
        }
        update_reward_full(i);
        update_predicted_reward(i);
        update_DA(i);
        update_w_ctx_msn(i);
        update_w_pf_tan(i);
        record_trial(jj,i);
        trial++;
    }
}

void simulate_acquisition_partial(int simulation) {
    int jj = simulation;
    for(int i=0; i<num_acquisition_trials; i++) {
        if(trial % 10 == 0) printf("Simulating trial %i\n", trial);
        reset_trial();
        for(int j=1; j<num_steps; j++) {
            update_pf(j);
            update_tan(j);
            update_sensory(j);
            update_msn(j);
            update_motor(j);
            update_response(j);
//            if(response) break;
        }
        update_reward_partial(i);
        update_predicted_reward(i);
        update_DA(i);
        update_w_ctx_msn(i);
        update_w_pf_tan(i);
        record_trial(jj,i);
        trial++;
    }
}

void simulate_extinction(int simulation) {
    int jj = simulation;
    for(int i=0; i<num_extinction_trials; i++) {
        if(trial % 10 == 0) printf("Simulating trial %i\n", trial);
        reset_trial();
        for(int j=1; j<num_steps; j++) {
            update_pf(j);
            update_tan(j);
            update_sensory(j);
            update_msn(j);
            update_motor(j);
            update_response(j);
//            if(response) break;
        }
        update_reward_extinction(i);
        update_predicted_reward(i);
        update_DA(i);
        update_w_ctx_msn(i);
        update_w_pf_tan(i);
        record_trial(jj,i);
        trial++;
    }
}

void initialize_parameters(void) {
    tau = 1.0;
    
    num_steps = 3000;
    num_acquisition_trials = 300;
    num_extinction_trials = 300;
    num_trials = num_acquisition_trials + num_extinction_trials;
    num_simulations = 100;
    
    trial=0;
    simulation=0;
    
    cue_onset = 1000;
    cue_duration = 1000;
    pf_cue_duration = cue_duration;
    
    srand((int)time(NULL));
    
    alpha_func_a = 1.0;
    alpha_func_b = 100.0;
    spike_length = floor(7.64*alpha_func_b);
    
    alpha_func_a_camkII = 50.0;
    alpha_func_b_camkII = 150.0;
    spike_length_camkII = floor(7.64*alpha_func_b_camkII);
    
    sensory_amp = 1500.0;
    pf_amp = 350.0;
    pause_mod_amp = 2.7;
    pause_decay = 0.0018;
    
    w_tan_msn = 125.0;
    w_msn_mot = 1.0;
    
    w_pf_tan_init = 0.2;
    w_ctx_msn_init = 0.2;
    
    pr_alpha = 0.01;
    
    response_threshold = 15.0;
    AMPA_threshold = 1.0;
    NMDA_threshold = 10.0;
    DA_base = 0.2;
    LTP_msn = 0.4e-6;
    LTD_msn = 0.5e-7;
    LTP_tan = 0.9e-6;
    LTD_tan = 0.1e-6;
}

void allocate_buffers(void) {
    spike = (float*) calloc(spike_length, sizeof(float));
    spike_camkII = (float*) calloc(spike_length, sizeof(float));
    msn_camkII = (float *) calloc(num_steps, sizeof(float));
    tan_camkII = (float *) calloc(num_steps, sizeof(float));
    
    sensory_square = (float *) calloc(num_steps, sizeof(float));
    
    msn_v = (float *) calloc(num_steps, sizeof(float));
    msn_u = (float *) calloc(num_steps, sizeof(float));
    msn_output = (float *) calloc(num_steps, sizeof(float));
    msn_spikes = (int *) calloc(num_steps, sizeof(int));
    
    pause_mod = (float *) calloc(num_steps, sizeof(float));
    pf_square = (float *) calloc(num_steps, sizeof(float));
    
    tan_v = (float *) calloc(num_steps, sizeof(float));
    tan_u = (float *) calloc(num_steps, sizeof(float));
    tan_output = (float *) calloc(num_steps, sizeof(float));
    tan_spikes = (int *) calloc(num_steps, sizeof(int));
    
    motor_v = (float *) calloc(num_steps, sizeof(float));
    motor_output = (float *) calloc(num_steps, sizeof(float));
    motor_spikes = (int *) calloc(num_steps, sizeof(int));
    
    // Records
    w_ctx_msn_record = (float **) calloc(num_simulations, sizeof(float*));
    w_pf_tan_record = (float **) calloc(num_simulations, sizeof(float*));
    response_record = (int **) calloc(num_simulations, sizeof(int*));
    response_time_record = (int **) calloc(num_simulations, sizeof(int*));
    reward_record = (int **) calloc(num_simulations, sizeof(int*));
    predicted_reward_record = (float **) calloc(num_simulations, sizeof(float*));
    rpe_record = (float **) calloc(num_simulations, sizeof(float*));
    DA_record = (float **) calloc(num_simulations, sizeof(float*));
    
    sensory_square_record = (float ***) calloc(num_simulations, sizeof(float**));
    msn_v_record = (float ***) calloc(num_simulations, sizeof(float**));
    msn_output_record = (float ***) calloc(num_simulations, sizeof(float**));
    msn_spikes_record = (int ***) calloc(num_simulations, sizeof(int**));
    pf_square_record = (float ***) calloc(num_simulations, sizeof(float**));
    tan_v_record = (float ***) calloc(num_simulations, sizeof(float**));
    tan_output_record = (float ***) calloc(num_simulations, sizeof(float**));
    tan_spikes_record = (int ***) calloc(num_simulations, sizeof(int**));
    motor_v_record = (float ***) calloc(num_simulations, sizeof(float**));
    motor_output_record = (float ***) calloc(num_simulations, sizeof(float**));
    motor_spikes_record = (int ***) calloc(num_simulations, sizeof(int**));
    
    for(int i=0; i<num_simulations; i++) {
        w_ctx_msn_record[i] = (float *) calloc(num_trials, sizeof(float));
        w_pf_tan_record[i] = (float *) calloc(num_trials, sizeof(float));
        response_record[i] = (int *) calloc(num_trials, sizeof(int));
        response_time_record[i] = (int *) calloc(num_trials, sizeof(float));
        reward_record[i] = (int *) calloc(num_trials, sizeof(int));
        predicted_reward_record[i] = (float *) calloc(num_trials, sizeof(float));
        rpe_record[i] = (float *) calloc(num_trials, sizeof(float));
        DA_record[i] = (float *) calloc(num_trials, sizeof(float));
        
        sensory_square_record[i] = (float **) calloc(num_trials, sizeof(float*));
        msn_v_record[i] = (float **) calloc(num_trials, sizeof(float*));
        msn_output_record[i] = (float **) calloc(num_trials, sizeof(float*));
        msn_spikes_record[i] = (int **) calloc(num_trials, sizeof(int*));
        pf_square_record[i] = (float **) calloc(num_trials, sizeof(float*));
        tan_v_record[i] = (float **) calloc(num_trials, sizeof(float*));
        tan_output_record[i] = (float **) calloc(num_trials, sizeof(float*));
        tan_spikes_record[i] = (int **) calloc(num_trials, sizeof(int*));
        motor_v_record[i] = (float **) calloc(num_trials, sizeof(float*));
        motor_output_record[i] = (float **) calloc(num_trials, sizeof(float*));
        motor_spikes_record[i] = (int **) calloc(num_trials, sizeof(int*));
    }
    
    for(int i=0; i<num_simulations; i++) {
        for(int j=0; j<num_trials; j++) {
            sensory_square_record[i][j] = (float *) calloc(num_steps, sizeof(float));
            msn_v_record[i][j] = (float *) calloc(num_steps, sizeof(float));
            msn_output_record[i][j] = (float *) calloc(num_steps, sizeof(float));
            msn_spikes_record[i][j] = (int *) calloc(num_steps, sizeof(int));
            pf_square_record[i][j] = (float *) calloc(num_steps, sizeof(float));
            tan_v_record[i][j] = (float *) calloc(num_steps, sizeof(float));
            tan_output_record[i][j] = (float *) calloc(num_steps, sizeof(float));
            tan_spikes_record[i][j] = (int *) calloc(num_steps, sizeof(int));
            motor_v_record[i][j] = (float *) calloc(num_steps, sizeof(float));
            motor_output_record[i][j] = (float *) calloc(num_steps, sizeof(float));
            motor_spikes_record[i][j] = (int *) calloc(num_steps, sizeof(int));
        }
    }
    
    for (int i=0; i<num_steps; i++) {
        sensory_square_record[simulation][trial][i] = sensory_square[i];
        msn_v_record[simulation][trial][i] = msn_v[i];
        msn_output_record[simulation][trial][i] = msn_output[i];
        msn_spikes_record[simulation][trial][i] = msn_spikes[i];
        pf_square_record[simulation][trial][i] = pause_mod[i];
        tan_v_record[simulation][trial][i] = tan_v[i];
        tan_output_record[simulation][trial][i] = tan_output[i];
        tan_spikes_record[simulation][trial][i] = tan_spikes[i];
        motor_v_record[simulation][trial][i] = motor_v[i];
        motor_output_record[simulation][trial][i] = motor_output[i];
        motor_spikes_record[simulation][trial][i] = motor_spikes[i];
    }
    
    // Average records
    w_ctx_msn_record_ave = (float *) calloc(num_trials, sizeof(float));
    w_pf_tan_record_ave = (float *) calloc(num_trials, sizeof(float));
    response_record_ave = (float *) calloc(num_trials, sizeof(float));
    response_time_record_ave = (float *) calloc(num_trials, sizeof(float));
    reward_record_ave = (float *) calloc(num_trials, sizeof(float));
    predicted_reward_record_ave = (float *) calloc(num_trials, sizeof(float));
    rpe_record_ave = (float *) calloc(num_trials, sizeof(float));
    DA_record_ave = (float *) calloc(num_trials, sizeof(float));
    
    sensory_square_record_ave = (float **) calloc(num_trials, sizeof(float*));
    msn_v_record_ave = (float **) calloc(num_trials, sizeof(float*));
    msn_output_record_ave = (float **) calloc(num_trials, sizeof(float*));
    msn_spikes_record_ave = (float **) calloc(num_trials, sizeof(float*));
    pf_square_record_ave = (float **) calloc(num_trials, sizeof(float*));
    tan_v_record_ave = (float **) calloc(num_trials, sizeof(float*));
    tan_output_record_ave = (float **) calloc(num_trials, sizeof(float*));
    tan_spikes_record_ave = (float **) calloc(num_trials, sizeof(float*));
    motor_v_record_ave = (float **) calloc(num_trials, sizeof(float*));
    motor_output_record_ave = (float **) calloc(num_trials, sizeof(float*));
    motor_spikes_record_ave = (float **) calloc(num_trials, sizeof(float*));
    
    for(int i=0; i<num_trials; i++) {
        sensory_square_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        msn_v_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        msn_output_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        msn_spikes_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        pf_square_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        tan_v_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        tan_output_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        tan_spikes_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        motor_v_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        motor_output_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
        motor_spikes_record_ave[i] = (float *) calloc(num_steps, sizeof(float));
    }
}

void initialize_buffers(void) {
    // Initialize spike output
    for(int i=0; i<spike_length; i++) {
		spike[i] = alpha_func_a*((float)i/alpha_func_b)*exp(-1*(i-alpha_func_b)/alpha_func_b);
	}
    
    for(int i=0; i<spike_length_camkII; i++) {
		spike_camkII[i] = alpha_func_a_camkII*((float)i/alpha_func_b_camkII)*exp(-1*(i-alpha_func_b_camkII)/alpha_func_b_camkII);
	}
}

void record_trial(int sim, int trl) {
    w_ctx_msn_record[simulation][trial] = w_ctx_msn;
    w_pf_tan_record[simulation][trial] = w_pf_tan;
    response_record[simulation][trial] = response;
    response_time_record[simulation][trial] = response_time;
    
    reward_record[simulation][trial] = reward;
    predicted_reward_record[simulation][trial] = predicted_reward;
    rpe_record[simulation][trial] = rpe;
    DA_record[simulation][trial] = DA;
    
    for (int i=0; i<num_steps; i++) {
        sensory_square_record[simulation][trial][i] = sensory_square[i];
        msn_v_record[simulation][trial][i] = msn_v[i];
        msn_output_record[simulation][trial][i] = msn_output[i];
        msn_spikes_record[simulation][trial][i] = msn_spikes[i];
        pf_square_record[simulation][trial][i] = pause_mod[i];
        tan_v_record[simulation][trial][i] = tan_v[i];
        tan_output_record[simulation][trial][i] = tan_output[i];
        tan_spikes_record[simulation][trial][i] = tan_spikes[i];
        motor_v_record[simulation][trial][i] = motor_v[i];
        motor_output_record[simulation][trial][i] = motor_output[i];
        motor_spikes_record[simulation][trial][i] = motor_spikes[i];
    }

//    memcpy(sensory_square_record[simulation][trial], sensory_square, num_steps*sizeof(float));
//    memcpy(msn_v_record[simulation][trial], msn_v, num_steps*sizeof(float));
//    memcpy(msn_output_record[simulation][trial], msn_output, num_steps*sizeof(float));
//    memcpy(msn_spikes_record[simulation][trial], msn_spikes, num_steps*sizeof(int));
//    memcpy(pf_square_record[simulation][trial], pause_mod, num_steps*sizeof(float));
//    memcpy(tan_v_record[simulation][trial], tan_v, num_steps*sizeof(float));
//    memcpy(tan_output_record[simulation][trial], tan_output, num_steps*sizeof(float));
//    memcpy(tan_spikes_record[simulation][trial], tan_spikes, num_steps*sizeof(int));
//    memcpy(motor_v_record[simulation][trial], motor_v, num_steps*sizeof(float));
//    memcpy(motor_output_record[simulation][trial], motor_output, num_steps*sizeof(float));
//    memcpy(motor_spikes_record[simulation][trial], motor_spikes, num_steps*sizeof(int));
}

void reset_trial(void) {
    response = 0;
    response_time = 0;
    reward = 0;
    rpe = 0;
    DA = 0;
    // NOTE: do not reset predicted_reward
    
    memset(msn_camkII, 0, num_steps*sizeof(float));
    memset(tan_camkII, 0, num_steps*sizeof(float));
    memset(sensory_square, 0, num_steps*sizeof(float));
    memset(msn_v, 0, num_steps*sizeof(float));
    memset(msn_u, 0, num_steps*sizeof(float));
    memset(msn_output, 0, num_steps*sizeof(float));
    memset(msn_spikes, 0, num_steps*sizeof(int));
    memset(pf_square, 0, num_steps*sizeof(float));
    memset(pf_square, 0, num_steps*sizeof(float));
    memset(pause_mod, 0, num_steps*sizeof(float));
    memset(tan_v, 0, num_steps*sizeof(float));
    memset(tan_u, 0, num_steps*sizeof(float));
    memset(tan_output, 0, num_steps*sizeof(float));
    memset(tan_spikes, 0, num_steps*sizeof(int));
    memset(motor_v, 0, num_steps*sizeof(float));
    memset(motor_output, 0, num_steps*sizeof(float));
    memset(motor_spikes, 0, num_steps*sizeof(int));
}

void reset_simulation(void) {
    trial = 0;
    w_pf_tan = w_pf_tan_init;
    w_ctx_msn = w_ctx_msn_init;
}

void reset_records(void) {
    memset(w_ctx_msn_record_ave, 0, num_trials*sizeof(float));
    memset(w_pf_tan_record_ave, 0, num_trials*sizeof(float));
    memset(response_record_ave, 0, num_trials*sizeof(float));
    memset(response_time_record_ave, 0, num_trials*sizeof(float));
    memset(reward_record_ave, 0, num_trials*sizeof(float));
    memset(predicted_reward_record_ave, 0, num_trials*sizeof(float));
    memset(rpe_record_ave, 0, num_trials*sizeof(float));
    memset(DA_record_ave, 0, num_trials*sizeof(float));
    
    for(int i=0; i<num_simulations; i++) {
        memset(w_ctx_msn_record[i], 0, num_trials*sizeof(float));
        memset(w_pf_tan_record[i], 0, num_trials*sizeof(float));
        memset(response_record[i], 0, num_trials*sizeof(int));
        memset(response_time_record[i], 0, num_trials*sizeof(int));
        memset(reward_record[i], 0, num_trials*sizeof(int));
        memset(predicted_reward_record[i], 0, num_trials*sizeof(float));
        memset(rpe_record[i], 0, num_trials*sizeof(float));
        memset(DA_record[i], 0, num_trials*sizeof(float));
        
        for (int j=0; j<num_trials; j++) {
            memset(sensory_square_record[i][j], 0, num_steps*sizeof(float));
            memset(msn_v_record[i][j], 0, num_steps*sizeof(float));
            memset(msn_output_record[i][j], 0, num_steps*sizeof(float));
            memset(msn_spikes_record[i][j], 0, num_steps*sizeof(int));
            memset(pf_square_record[i][j], 0, num_steps*sizeof(float));
            memset(tan_v_record[i][j], 0, num_steps*sizeof(float));
            memset(tan_output_record[i][j], 0, num_steps*sizeof(float));
            memset(tan_spikes_record[i][j], 0, num_steps*sizeof(int));
            memset(motor_v_record[i][j], 0, num_steps*sizeof(float));
            memset(motor_output_record[i][j], 0, num_steps*sizeof(float));
            memset(motor_spikes_record[i][j], 0, num_steps*sizeof(int));
        }
    }
}

void compute_average_results(void) {
    for(int i=0; i<num_simulations; i++) {
        for(int j=0; j<num_trials; j++) {
            w_ctx_msn_record_ave[j] += w_ctx_msn_record[i][j];
            w_pf_tan_record_ave[j] += w_pf_tan_record[i][j];
            response_record_ave[j] += response_record[i][j];
            response_time_record_ave[j] += response_time_record[i][j];
            reward_record_ave[j] += reward_record[i][j];
            predicted_reward_record_ave[j] += predicted_reward_record[i][j];
            rpe_record_ave[j] += rpe_record[i][j];
            DA_record_ave[j] += DA_record[i][j];
            
            for(int k=0; k<num_steps; k++) {
                sensory_square_record_ave[j][k] += sensory_square_record[i][j][k];
                msn_v_record_ave[j][k] += msn_v_record[i][j][k];
                msn_output_record_ave[j][k] += msn_output_record[i][j][k];
                msn_spikes_record_ave[j][k] += msn_spikes_record[i][j][k];
                pf_square_record_ave[j][k] += pf_square_record[i][j][k];
                tan_v_record_ave[j][k] += tan_v_record[i][j][k];
                tan_output_record_ave[j][k] += tan_output_record[i][j][k];
                tan_spikes_record_ave[j][k] += tan_spikes_record[i][j][k];
                motor_v_record_ave[j][k] += motor_v_record[i][j][k];
                motor_output_record_ave[j][k] += motor_output_record[i][j][k];
                motor_spikes_record_ave[j][k] += motor_spikes_record[i][j][k];
            }
        }
    }
    
    for(int j=0; j<num_trials; j++) {
        w_ctx_msn_record_ave[j] /= (float) num_simulations;
        w_pf_tan_record_ave[j] /= (float) num_simulations;
        response_record_ave[j] /= (float) num_simulations;
        response_time_record_ave[j] /= (float) num_simulations;
        
        reward_record_ave[j] /= (float) num_simulations;
        predicted_reward_record_ave[j] /= (float) num_simulations;
        rpe_record_ave[j] /= (float) num_simulations;
        DA_record_ave[j] /= (float) num_simulations;
        
        for(int k=0; k<num_steps; k++) {
            sensory_square_record_ave[j][k] /= (float) num_simulations;
            msn_v_record_ave[j][k] /= (float) num_simulations;
            msn_output_record_ave[j][k] /= (float) num_simulations;
            msn_spikes_record_ave[j][k] /= (float) num_simulations;
            pf_square_record_ave[j][k] /= (float) num_simulations;
            tan_v_record_ave[j][k] /= (float) num_simulations;
            tan_output_record_ave[j][k] /= (float) num_simulations;
            tan_spikes_record_ave[j][k] /= (float) num_simulations;
            motor_v_record_ave[j][k] /= (float) num_simulations;
            motor_output_record_ave[j][k] /= (float) num_simulations;
            motor_spikes_record_ave[j][k] /= (float) num_simulations;
        }
    }
}

void save_average_essentials(void) {
    FILE *fid_sensory_output;
    FILE *fid_msn_output;
    FILE *fid_pf_output;
    FILE *fid_tan_output;
    FILE *fid_motor_output;
    FILE *fid_response;
    FILE *fid_response_time;
    FILE *fid_w_ctx_msn;
    FILE *fid_w_pf_tan;
    FILE *fid_reward;
    FILE *fid_predicted_reward;
    FILE *fid_rpe;
    FILE *fid_DA;
    
//    fid_sensory_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_sensory_output_0", "w");
//    fid_msn_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_msn_output_0", "w");
//    fid_pf_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_pf_output_0", "w");
//    fid_tan_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_tan_output_0", "w");
//    fid_motor_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_motor_output_0", "w");
//    fid_response = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_response_0", "w");
//    fid_response_time = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_response_time_0", "w");
//    fid_w_ctx_msn = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_w_ctx_msn_0", "w");
//    fid_w_pf_tan = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_w_pf_tan_0", "w");
//    fid_reward = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_reward_0", "w");
//    fid_predicted_reward = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_predicted_reward_0", "w");
//    fid_rpe = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_rpe_0", "w");
//    fid_DA = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE_DA_0", "w");
    
    switch (condition_flag) {
            
        case 0:
            fid_sensory_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/sensory_output_0", "w");
            fid_msn_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/msn_output_0", "w");
            fid_pf_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/pf_output_0", "w");
            fid_tan_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/tan_output_0", "w");
            fid_motor_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/motor_output_0", "w");
            fid_response = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/response_0", "w");
            fid_response_time = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/response_time_0", "w");
            fid_w_ctx_msn = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/w_ctx_msn_0", "w");
            fid_w_pf_tan = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/w_pf_tan_0", "w");
            fid_reward = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/reward_0", "w");
            fid_predicted_reward = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/predicted_reward_0", "w");
            fid_rpe = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/rpe_0", "w");
            fid_DA = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/DA_0", "w");
            break;
            
        case 1:
            fid_sensory_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/sensory_output_1", "w");
            fid_msn_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/msn_output_1", "w");
            fid_pf_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/pf_output_1", "w");
            fid_tan_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/tan_output_1", "w");
            fid_motor_output = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/motor_output_1", "w");
            fid_response = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/response_1", "w");
            fid_response_time = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/response_time_1", "w");
            fid_w_ctx_msn = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/w_ctx_msn_1", "w");
            fid_w_pf_tan = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/w_pf_tan_1", "w");
            fid_reward = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/reward_1", "w");
            fid_predicted_reward = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/predicted_reward_1", "w");
            fid_rpe = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/rpe_1", "w");
            fid_DA = fopen("/Users/mjc/Documents/mjc/instrumental_discrete/instrumental_discrete_PRE/DA_1", "w");
            break;
            
        default:
            break;
    }
    
    for(int i=0; i<num_trials; i++) {
        for(int j=0; j<num_steps; j++) {
            fprintf(fid_sensory_output, "%g\n", sensory_square_record_ave[i][j]);
            fprintf(fid_msn_output, "%g\n",msn_output_record_ave[i][j]);
            fprintf(fid_pf_output, "%g\n", pf_square_record_ave[i][j]);
            fprintf(fid_tan_output, "%g\n",tan_output_record_ave[i][j]);
            fprintf(fid_motor_output, "%g\n",motor_output_record_ave[i][j]);
        }
        
        fprintf(fid_w_ctx_msn, "%g\n",w_ctx_msn_record_ave[i]);
        fprintf(fid_w_pf_tan, "%g\n",w_pf_tan_record_ave[i]);
        fprintf(fid_response, "%g\n",response_record_ave[i]);
        fprintf(fid_response_time, "%g\n",response_time_record_ave[i]);
        fprintf(fid_reward, "%g\n",reward_record_ave[i]);
        fprintf(fid_predicted_reward, "%g\n",predicted_reward_record_ave[i]);
        fprintf(fid_rpe, "%g\n",rpe_record_ave[i]);
        fprintf(fid_DA, "%g\n",DA_record_ave[i]);
    }
    
    fclose(fid_sensory_output);
    fclose(fid_msn_output);
    fclose(fid_pf_output);
    fclose(fid_tan_output);
    fclose(fid_motor_output);
    fclose(fid_response);
    fclose(fid_response_time);
    fclose(fid_w_ctx_msn);
    fclose(fid_w_pf_tan);
    fclose(fid_reward);
    fclose(fid_predicted_reward);
    fclose(fid_rpe);
    fclose(fid_DA);
}

void update_pf(int tick) {
    int t = tick;
    if(t >= cue_onset && t < cue_onset+pf_cue_duration) {
        pf_square[t] = pf_amp;
        pause_mod[t] = pf_amp;
    }
    
    if(t >= cue_onset+pf_cue_duration) {
        pause_mod[t] = pf_amp*exp(-pause_decay*(t-(cue_onset+pf_cue_duration)));
    }
}

void update_tan(int tick) {
    int t = tick;
    
	// TAN - modified intrinsically bursting
	float C=100.0, vr=-75.0, vt=-45.0, k=1.2;
	float a=0.01, b=5.0, c=-56.0, D=130.0;
	float vpeak=60.0;
	float E=950.0;
	
	tan_v[0] = vr;
    tan_v[1] = vr;
	
    noise = randn(0.0, 10.0);
    tan_v[t] = tan_v[t-1] + tau*(k*(tan_v[t-1]-vr)*(tan_v[t-1]-vt) - tan_u[t-1] + E + w_pf_tan*pf_square[t-1] + noise)/C;
    tan_u[t] = tan_u[t-1] + tau*a*(b*(tan_v[t-1]-vr) - tan_u[t-1] + w_pf_tan*pause_mod_amp*pause_mod[t-1]);
    
    if(tan_v[t]>=vpeak) {
        tan_v[t-1]=vpeak;
        tan_v[t]=c;
        tan_u[t]=tan_u[t]+D;
        
        // track cell output
        if(t<num_steps-spike_length) {
            for(int j=0; j<spike_length; j++) {
                tan_output[t+j] += spike[j];
            }
        } else {
            for(int j=0; j<num_steps-t; j++) {
                tan_output[t+j] += spike[j];
            }
        }
        
        if(t<num_steps-spike_length_camkII) {
            for(int j=0; j<spike_length_camkII; j++) {
                tan_camkII[t+j] += spike_camkII[j];
            }
        } else {
            for(int j=0; j<num_steps-t; j++) {
                tan_camkII[t+j] += spike_camkII[j];
            }
        }
    }
}

void update_sensory(int tick) {
    int t = tick;
    if(t >= cue_onset && t <= cue_onset+cue_duration) sensory_square[t-1] = sensory_amp;
}

void update_msn(int tick) {
    int t = tick;
    
	// MSN -- striatum
	float C=50, vr=-80, vt=-25, k=1;
	float a=0.01, b=-20, c=-55, D=150;
	float vpeak=40;
	
	msn_v[0] = vr;
    msn_v[1] = vr;
	
    noise = randn(125.0, 10.0);
    msn_v[t] = msn_v[t-1] + tau*(k*(msn_v[t-1]-vr)*(msn_v[t-1]-vt) - msn_u[t-1] + w_ctx_msn*pos(sensory_square[t-1]-w_tan_msn*tan_output[t-1]) + noise)/C;
    msn_u[t] = msn_u[t-1]+tau*a*(b*(msn_v[t-1]-vr)-msn_u[t-1]);
    
    if(msn_v[t]>=vpeak) {
        msn_spikes[t-1] = 1;
        msn_v[t-1]=vpeak;
        msn_v[t]=c;
        msn_u[t]=msn_u[t]+D;
        
        // track cell output
        if(t<num_steps-spike_length) {
            for(int j=0; j<spike_length; j++) {
                msn_output[t+j] += spike[j];
            }
        } else {
            for(int j=0; j<num_steps-t; j++) {
                msn_output[t+j] += spike[j];
            }
        }
        
        // track cell output
        if(t<num_steps-spike_length_camkII) {
            for(int j=0; j<spike_length_camkII; j++) {
                msn_camkII[t+j] += spike_camkII[j];
            }
        } else {
            for(int j=0; j<num_steps-t; j++) {
                msn_camkII[t+j] += spike_camkII[j];
            }
        }
    }
}

void update_motor(int tick) {
    int t = tick;
    
    // Regular Spiking Neuron
	float C=25, vr=-60, vt=-40, k=0.7;
	float c=-50;
	float vpeak=35;
	
	motor_v[0] = vr;
    motor_v[1] = vr;
    
    noise = randn(75.0, 50.0);
    motor_v[t] = motor_v[t-1] + tau*( k*(motor_v[t-1]-vr)*(motor_v[t-1]-vt) + w_msn_mot*msn_output[t-1] + noise)/C;
    
    if(motor_v[t]>=vpeak) {
        motor_spikes[t-1] = 1;
        motor_v[t]=vpeak;
        motor_v[t]=c;
        
        if(t<num_steps-spike_length) {
            for(int j=0; j<spike_length; j++)
            {
                motor_output[t+j] += spike[j];
            }
        }else {
            for(int j=0; j<num_steps-t; j++) {
                motor_output[t+j] += spike[j];
            }
        }
    }
}

void update_response(int tick) {
    if( motor_output[tick] >= response_threshold && response == 0 ) {
        response = 1;
        response_time = tick;
    }
    
    if(response == 0 && rand() / (float) RAND_MAX > 0.9999) {
        response = 1;
        response_time = tick;
    }
}

void update_reward_full(int trial) {
    if(response) reward = 1;
}

void update_reward_partial(int trial) {
    if(response && rand()/ (float) RAND_MAX < 0.5) reward = 1;
}

void update_reward_extinction(int trial) {
    reward = 0;
}

void update_predicted_reward(int trial) {
    rpe = reward - predicted_reward;
    predicted_reward += pr_alpha*rpe;
}

void update_DA(int trial) {
    DA = cap(0.8*rpe + DA_base);
}

void update_w_ctx_msn(int trial) {
    // Compute weight change
    w_delta = 0.0;
    for(int i=cue_onset; i<cue_onset+cue_duration; i++) {
        w_delta += LTP_msn*pos( msn_camkII[i] - NMDA_threshold )*pos( DA - DA_base )*cap( 1.0 - w_ctx_msn )
                 - LTD_msn*pos( msn_camkII[i] - NMDA_threshold )*pos( DA_base - DA )*cap( w_ctx_msn );
    }
    // Update w_ctx_msn
    w_ctx_msn += w_delta;
    w_ctx_msn = cap(w_ctx_msn);
}

void update_w_pf_tan(int trial) {
    // Compute weight change
    w_delta = 0.0;
    for(int i=cue_onset; i<cue_onset+cue_duration; i++) {
        w_delta += LTP_tan*pos( tan_camkII[i] - NMDA_threshold )*pos( DA - DA_base )*cap( 1.0 - w_pf_tan )
                 - LTD_tan*pos( tan_camkII[i] - NMDA_threshold )*pos( DA_base - DA )*cap( w_pf_tan );
    }
    // Update w_ctx_msn
    w_pf_tan += w_delta;
    w_pf_tan = cap(w_pf_tan);
}

# pragma mark -- Utilities --
float cap(float val) {
    if(val < 0.0) {
        val = 0.0;
    }else if(val > 1.0) {
        val = 1.0;
    }
    return val;
}

float pos(float val) {
    return val > 0.0 ? val : 0.0;
}

float randn(float mean, float variance) {
    float mu = mean;
    float sigma = sqrtf(variance);
    float uni_noise = (float) rand()/RAND_MAX;
    float normal_noise = mu - (sigma*sqrt(3.0)/3.14159)*(logf(1.0-uni_noise) - logf(uni_noise));
    return normal_noise;
}
