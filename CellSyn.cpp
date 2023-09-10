#include "CellSyn.h"


CellSyn::CellSyn(int i,int j,int b_size,int t_size,Cell_Info *cell_info,int type){    

  rand_buffer = new drand48_data;
  srand48_r(2000,rand_buffer);
  //printf("srand 48 \n");
  //srand48_r(time(NULL),rand_buffer);

  this->old_v=0;

  this->type = type;
  this->cell_info = cell_info;
  this->num_syns = this->cell_info->total_connects;

  num_spikes = 0;
    
  flip = 0;
  flip_for_maps = 0;

  m = i;
  n = j;
  lhs_size = t_size;
  y = new double[lhs_size];
  f = new double[lhs_size];
  s = new double[lhs_size];
  yk = new double[lhs_size];
  for(int iter = 0; iter < lhs_size; iter++){
	y[iter] = 0;
	f[iter] = 0;
	s[iter] = 0;
	yk[iter] = 0;
  }

  current_step = 0;
  cur_time = 0;
  xk = 0;
  this->tau = TAU;
  num_b_eq = b_size;
}

//TODO this function should use enum Cell_Type instead of plain numbers 
//this function creates all the cells(CellSyns) and does some basic initialization
//all arguments are inputs.
CellSyn* CellSyn::initialize(int type,int i, int j, Cell_Info ****cells_info,string output_location ){
  
  CellSyn *me;

  switch(type){    
  case E_RE:
    me = new REsyn(i,j,N_RE,N_RE,cells_info[E_RE][i][j],type);
    return me;

  case E_REa:
    me = new REsyn(i,j,N_RE,N_RE,cells_info[E_REa][i][j],type);
    return me;

  case E_TC:
    me = new TCcore(i,j,N_TC,N_TC+N_GB*cells_info[E_TC][i][j]->num_connects_norm[E_RE],cells_info[E_TC][i][j],type);
    return me;

  case E_TCa: 
    me = new TCmatrix(i,j,N_TCa,N_TCa+N_GB*cells_info[E_TCa][i][j]->num_connects_norm[E_REa],cells_info[E_TCa][i][j],type);
    return me;

  case E_CX:
	if (I_HH==1){
	  me = new CXsyn(i,j,N_CX,N_CX+N_GB*cells_info[E_CX][i][j]->num_connects_norm[E_IN],cells_info[E_CX][i][j],type);}
	else {
	  // printf("\n Creating maps");
	  me = new CXsyn_Map(i,j,N_CX,N_CX+N_GB*cells_info[E_CX][i][j]->num_connects_norm[E_IN],cells_info[E_CX][i][j],type);}
    return me;

  case E_CXa:
	if (I_HH==1){
	  me = new CXasyn(i,j,N_CXa,N_CXa+N_GB*cells_info[E_CXa][i][j]->num_connects_norm[E_INa],cells_info[E_CXa][i][j],type);}
	else{
	  me = new CXasyn_Map(i,j,N_CXa,N_CXa+N_GB*cells_info[E_CXa][i][j]->num_connects_norm[E_INa],cells_info[E_CXa][i][j],type);}
    return me;

  case E_IN:
	if (I_HH==1){
	  me = new INsynCore(i,j,N_IN,N_IN,cells_info[E_IN][i][j],type);}
	else{
	  me = new INsynCore_Map(i,j,N_IN,N_IN,cells_info[E_IN][i][j],type);}
    return me;

  case E_INa:
	if (I_HH==1){
	  me = new INsynMatrix(i,j,N_INa,N_INa,cells_info[E_INa][i][j],type);}
	else{
	  me = new INsynMatrix_Map(i,j,N_INa,N_INa,cells_info[E_INa][i][j],type); }
    return me;

  case E_CX6:
	if (I_HH==1){
	  me = new CXsyn6(i,j,N_CX6,N_CX6+N_GB*cells_info[E_CX6][i][j]->num_connects_norm[E_IN6],cells_info[E_CX6][i][j],type);}
	else{
	  me = new CXsyn6_Map(i,j,N_CX6,N_CX6+N_GB*cells_info[E_CX6][i][j]->num_connects_norm[E_IN6],cells_info[E_CX6][i][j],type);}
    return me;

  case E_IN6:
	if (I_HH==1){
	  me = new INsyn6(i,j,N_IN6,N_IN6,cells_info[E_IN6][i][j],type);}
	else{
	  me = new INsyn6_Map(i,j,N_IN6,N_IN6,cells_info[E_IN6][i][j],type);}
    return me;
  default:
    printf("tried to create a non existant cell type in initialize\n");
    exit(1);
  }
}


// Called from memory function - syns[i] = initialize_synapses(syn_info[i]);
Syn* CellSyn::initialize_synapses(Syn_Info *syn_info){

  Syn *syni=NULL;
	  
  if(syn_info->type == E_AMPA){
	syni = new AMPA(syn_info);
  }else if(syn_info->type == E_GiantAMPA){
  	syni = new Syn(syn_info);
	syni->giant = get_GiantSyn(syn_info->from_type,syn_info->type);
  }else if(syn_info->type == E_GiantAMPAMap){
  	syni = new Syn(syn_info);
	syni->giant = get_GiantSyn(syn_info->from_type,syn_info->type);
  }else if(syn_info->type == E_GABA_A){
	syni = new GABA_A(syn_info);
  }else if(syn_info->type == E_AMPA_D3){ 
	syni = new AMPA_D3(syn_info);
  }else if(syn_info->type == E_AMPA_D2){
	syni = new AMPA_D2(syn_info); 	  	
  }else if(syn_info->type == E_AMPA_D1){
	syni = new AMPA_D1(syn_info); 
  }else if(syn_info->type == E_NMDA_D1){
	syni = new NMDA_D1(syn_info); 
  }else if(syn_info->type == E_GABA_B){
	syni = new GB(syn_info); 
  }else if(syn_info->type == E_GABA_A_D2){
	syni = new GABA_A_D2(syn_info); 
  }else if(syn_info->type == E_GABAMap_A_D2){
	syni = new GABAAmapD1(syn_info); 
  }else if(syn_info->type == E_AMPAMap_D1){
	syni = new AMPAmapD1(syn_info);
  }else if(syn_info->type == E_NMDAMap_D1){
	syni = new NMDAmapD1(syn_info); 
  }else if(syn_info->type == E_GAP){
	syni = new GAP(syn_info);
  }else{
	printf("Synapse initialization failed: type %d from neuron: %d to neuron %d \n", syn_info->type, syn_info->from_type, syn_info->to_type);
	exit(1);
  }

  return syni;
}


GiantSyn* CellSyn::define_GiantSyn(enum Cell_Type cell, enum Syn_Type syn){
  GiantSyn *S=NULL;
  //TODO add case statement for more synapses
  if (syn == E_GiantAMPAMap)
	S = new GiantSynAMPAMap();
  else if (syn == E_GiantAMPA)
	S = new GiantSynAMPA();

  if (!S) {printf("Error Giant synapse initialization!\n");exit(1);}
  return S;
}


void CellSyn::signal_spike(){
    
  double message = base_cell->get_v_soma();

  if(this->old_v < 0 && message >= 0){
	this->old_v = message;
	flip = 1;
	flip_for_maps=1;
	num_spikes = num_spikes +1;
	print_spike_time(cur_time,this->type,this->m,this->n);
  }else{
	flip = 0;
	this->old_v = message;
	return;
  }
 
}

void CellSyn::reset_flip_maps(){
    
  if(flip_for_maps == 1){
	flip_for_maps=0;
  }
 
}

void CellSyn::zero_temps(){
  // not sure if it is need 
  int i = 0;
  for(i =0; i<lhs_size; i++){
	s[i]=0;
	yk[i]=0;
  }
}

void CellSyn::reset_y(){
  int i = 0;
  for(i = 0; i < lhs_size; ++i){ 
	y[i] = y[i] + (tau/6.0)*(s[i] + f[i]);
  }
}

void CellSyn::step(){

  if (this->ismap){
	this->calc(cur_time,y,f,2);
	cur_time=cur_time+TAU_Map;
  } 
  else{
	int i = 0;
	if(current_step == 0){
	  this->calc(cur_time,y,f,1);
	  for(i = 0; i < lhs_size; ++i){ 
		s[i] = f[i]; 
		yk[i] = y[i] + (tau/2.0)*f[i]; 
	  }
	  xk = cur_time + tau/2.0;
	}
  
	if(current_step ==1){
	  this->calc(xk,yk,f,0);
	  for(i = 0; i < lhs_size; ++i) { 
		s[i] = s[i] + 2.0*f[i]; 
		yk[i] = y[i] + (tau/2.0)*f[i]; 
	  }
	}

	if(current_step == 2){
	  this->calc(xk,yk,f,0);
	  for(i = 0; i < lhs_size; ++i) { 
		s[i] = s[i] + 2.0*f[i]; 
		yk[i] = y[i] + tau*f[i]; 
	  }
	  xk = cur_time + tau;
	}

	if(current_step == 3){
	  this->calc(xk,yk,f,0);
	  cur_time = cur_time + tau;
	}
    
	current_step = current_step +1;
	//TODO kill this condition
	if(current_step > 3){
	  current_step = 0;
	}
  }

}


extern CellSyn **cells;
extern Cell_Info ****cells_info;
//this function runs thought each synapse for a cell and calculates its current before calculating the cell itself
//step == 0 <=> RK 2,3,4 step;  step == 1 <=> RK 1; step == 2 <=> MAP
void CellSyn::calc(double x, double *y_ini, double *f_ini, int step){
  
  double current = 0.0;
  int input_soma = 0; 	  // input_soma - pre-synaptic spike
  double input_dend = 0;  
  int i = 0;
  int k = 0;

  double vd=base_cell->get_v_dend();

  // reset giant synapse when step
  if (step>0)
    reset_GiantSyns();

  // Run through all synapses
  for(i =0; i < num_syns; i++){

	int syn_type=syns[i]->type;

	int from=syns[i]->from_cell;
	GiantSyn *gs = syns[i]->giant; // Check if this synapse part of giant synapse
	  
	if(syn_type == E_GAP){
	  input_dend = cells[from]->base_cell->get_v_dend();
	  current = current - ((GAP*)syns[i])->calc_gap(x,vd, input_dend); // /max_connect_gap[type][ftype];
	  continue; 
	}

	if (step>0) {
	  if (this->ismap)
		input_soma = cells[from]->flip_for_maps;
	  else {
		input_soma = cells[from]->flip;	// check if input neuron spiked
               }

	  // if it is giant synapse add spikes
	  if (gs && input_soma)
	    gs->add_spike();
    } else {
       input_soma=0; // Ensures the spike is used only when in step>0, since the reset occurs at step=1 or 2
    }

	if(gs){}
	else{
	  if(syn_type == E_GABA_B){
		//current = current - ((GB*)syns[i])->calc_gb(x,y+num_b_eq+2*k, f+num_b_eq+2*k, vd, input_soma)/max_connect[type][ftype];
		// scaling now implimented in initialization of network - Prevents division every time calc is called
		current = current - ((GB*)syns[i])->calc_gb(x,y+num_b_eq+2*k, f+num_b_eq+2*k, vd, input_soma);
		k = k + 1;
	  }
	  else if(syn_type==E_AMPA_D3)
	    current = current - ((AMPA_D3*)syns[i])->calc_ampad3(x,vd,input_soma,this->flip);
        else{ 
		//general case - for all other synapses
		// current = current - (syns[i])->calc(x, vd, input_soma)/max_connect[type][ftype];
		current = current - (syns[i])->calc(x, vd, input_soma);
	  }
	}
  }

  int period=1000;
  double stime=fmod(x,period);

  int period2=1333; //periodic stimulation
  double stime2=fmod(x,period2); 

  //  double stime2=fmod(x,2*period);

  // int gap=25;
  if((x>=10001&&x<=60000)||(x>=160001&&x<awake_end)||(x>stage3_end)){
    if ((stime>=0)&&(stime<10)) {
     if (this->cell_type == E_CX && this->m>=input_start_neuron+0 && this->m<=input_end_neuron+0){
        current=current + input_strength;
      }
     if (this->cell_type == E_CX && this->m>=input_start_neuron+96 && this->m<=input_end_neuron+96){
       // current=current + input_strength;
      }
    }
  }

  //awake seq1 training sessions
   if(x>=60001&&x<140000){
    if ((stime>=0)&&(stime<10)) {
      if (this->cell_type == E_CXa && this->m>=input_start_neuron && this->m<=input_end_neuron){
        current=current + input_strength;
      }
    }
    if ((stime>=15)&&(stime<25)) {
      if (this->cell_type == E_CX && this->m>=input_start_neuron+0 && this->m<=input_end_neuron+0){
        current=current + input_strength;
      }
     }
    if ((stime>=30)&&(stime<40)) {
      if (this->cell_type == E_CX && this->m>=input_start_neuron+5 && this->m<=input_end_neuron+5){
        current=current + input_strength;
      }
    }
   if ((stime>=45)&&(stime<55)) {
      if (this->cell_type == E_CX && this->m>=input_start_neuron+10 && this->m<=input_end_neuron+10){
        current=current + input_strength;
      }
    }

    if ((stime>=60)&&(stime<70)) {
      if (this->cell_type == E_CX && this->m>=input_start_neuron+15 && this->m<=input_end_neuron+15){
        current=current + input_strength;
      }
    }
    if ((stime>=75)&&(stime<85)) {
      if (this->cell_type == E_CX && this->m>=input_start_neuron+20 && this->m<=input_end_neuron+20){
        current=current + input_strength;
      }
    }
   }

   // awake seq2 training sessions
   /*  if(x>=60001&&x<140000){
     if ((stime>=0)&&(stime<10)) {
       if (this->cell_type == E_CXa && this->m>=input_start_neuron+96 && this->m<=input_end_neuron+96){
	 current=current + input_strength;
       }
     }
     if ((stime>=15)&&(stime<25)) {
       if (this->cell_type == E_CX && this->m>=input_start_neuron+96 && this->m<=input_end_neuron+96){
	 current=current + input_strength;
       }
     }
     if ((stime>=30)&&(stime<40)) {
       if (this->cell_type == E_CX && this->m>=input_start_neuron+91 && this->m<=input_end_neuron+91){
	 current=current + input_strength;
       }
     }
     if ((stime>=45)&&(stime<55)) {
       if (this->cell_type == E_CX && this->m>=input_start_neuron+86 && this->m<=input_end_neuron+86){
	 current=current + input_strength;
       }
     }

     if ((stime>=60)&&(stime<70)) {
       if (this->cell_type == E_CX && this->m>=input_start_neuron+81 && this->m<=input_end_neuron+81){
	 current=current + input_strength;
       }
     }
     if ((stime>=75)&&(stime<85)) {
       if (this->cell_type == E_CX && this->m>=input_start_neuron+76 && this->m<=input_end_neuron+76){
	 current=current + input_strength;
       }
     }
   }
   */


   // phase stimulation
   if(stimType==2) {   // stim during down state
     if(avgLFP<=-64 && abs(avgLFP+64)<2 &&lock==0){
	  lock=1;
	  stim_time0=x;
	  stim_time=x+dtime;
	}
      if(avgLFP>-60 && x>=stim_time0+400 && x>stim_time) {lock=0;}    
    }
   else if(stimType==3)         // stim during up state
    {
      if(avgLFP>=-64 && abs(avgLFP+64)<2 &&lock==0){
	  lock=1;
	  stim_time0=x;
	  stim_time=x+dtime;
	}
      if(avgLFP<-67 && x>=stim_time0+400 && x>stim_time) {lock=0;}    
    }else if(stimType==1)
     {
        // periodic stimulation
       if(stime2<0.01)
       {
          stim_time=x;
       }
     }
   else //No Stim
   {
       stim_time = -1 ;
    }

  // Cue presented during sleep
  if(x>=stim_start&&x<stim_end){
    if((x>stim_time)&&(x<stim_time+50)){
      if (this->cell_type == E_CXa && this->m>=input_start_neuron && this->m<=input_end_neuron){
	current=current + input_strength*0.3;//*0.3;
      }
     }
      if((x>stim_time)&&(x<stim_time+50)){
	if (this->cell_type == E_CXa && this->m>=input_start_neuron+96 && this->m<=input_end_neuron+96){
	  // current=current + input_strength*0.0;//*0.3;                                                                                     
	}
    }
  }

  current =  current - calc_GiantSyns();
  base_cell->calc(x, current, y_ini, f_ini);

}

GiantSyn* CellSyn::find_GiantSyn(enum Cell_Type cell, enum Syn_Type syn){
  list<GiantSyn*>::iterator it=giant_syns.begin();
  list<GiantSyn*>::iterator end=giant_syns.end();
  for (;it!=end;it++)
	if ( (*it)->from_type == cell && (*it)->syn_type == syn )
	  return (*it);
  return NULL;
}

GiantSyn* CellSyn::get_GiantSyn(enum Cell_Type cell, enum Syn_Type syn){
  //Here can be placed #ifndef GiantSyn return NULL; 
  GiantSyn *S=NULL;
  S = find_GiantSyn(cell,syn);
  if (!S) {
	S = define_GiantSyn(cell,syn);
	giant_syns.push_back(S);
  }
  return S;
}

void CellSyn::reset_GiantSyns(){
  list<GiantSyn*>::iterator it=giant_syns.begin();
  list<GiantSyn*>::iterator end=giant_syns.end();
  for (;it!=end;it++)
	(*it)->reset_spikes();
}

double CellSyn::calc_GiantSyns(){
  double sum = 0;
  list<GiantSyn*>::iterator it=giant_syns.begin();
  list<GiantSyn*>::iterator end=giant_syns.end();
  for (;it!=end;it++)
	sum += (*it)->calc();

  return sum;
}

