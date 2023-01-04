#
#
#
#

clear all; close all;

graphics_toolkit qt

global specs;

specs.ngrid = 100;
specs.viscosity = 1.0;
specs.timestep = 1.0e-5;
specs.numloops = 1000;
specs.samplefreq = 100;
specs.vortex1 = [-1.0, 0.0, 1.0];
specs.vortex2 = [1.0, 0.0, 1.0];
specs.vortex3 = [0.0, 0.0, 0.0];
specs.vortex4 = [0.0, 0.0, 0.0];
specs.surface = false;
specs.removefiles = false;
specs.euler = true;
specs.lbox = 10.0;



##
## Ad hoc boundary conditions
##
function phi = bc(phi)

	phi(1,:)   =0;
	phi(end,:) = 0;
	phi(:,1)   = 0;
	phi(:,end) = 0;

endfunction


##
## Simple explicit Euler -  direct method for the Poisson equation
##
function [w u v phi]  = eulerstep(w, u, v, cmatrix, invRe, dt, X, Y)
	
	[ngy, ngx] = size(X);
	
	[dxw dyw dx2w dy2w] = der2d(w, X, Y, 'n');

	f = -u.*dxw - v.*dyw + invRe.*(dx2w+dy2w);
	
	w = w + f*dt;
	w = bc(w);
	
	wvec  = -mattovec(w);
	phivec = cmatrix\wvec;
  phi  = vectomat(phivec, ngy, ngx);
	phi  = bc(phi);
		
	[v u]  = der2d(phi, X, Y, 'n');
	v      = -v;
	u = bc(u); v = bc(v);
		
endfunction
			  
function [f w u v]  = rkf(wn, w, un, vn, cmatrix, invRe, dt, X, Y, opt)
  
  [ngy, ngx] = size(X);

  [dxw dyw dx2w dy2w] = der2d(w, X, Y, 'n');
     
  f = -un.*dxw - vn.*dyw + invRe.*(dx2w+dy2w);

  if ( opt == 1 )
     w = wn + f*dt;
     w = bc(w);
     
     wvec   = -mattovec(w);

     phivec = cmatrix\wvec;
     phi = vectomat(phivec, ngy, ngx);
     phi = bc(phi);

     [v u] = der2d(phi, X, Y, 'n');
     v = -v;
     u = bc(u); v=bc(v);

  endif
  
endfunction

function [w, u, v, phi] =  rk4step(w, u, v, cmatrix, invRe, dt, X, Y)
  
  [ngy ngx] = size(X);
  
  wn = w;     

  ## f1  
  [f1 w1 u1 v1]  = rkf(wn, w, u, v, cmatrix, invRe, 0.5*dt, X, Y, 1);	

  ## f2
  [f2 w2 u2 v2]  = rkf(wn, w1, u1, v1, cmatrix, invRe, 0.5*dt, X, Y, 1);	
  
  ## f3
  [f3 w3 u3 v3]  = rkf(wn, w2, u2, v2, cmatrix, invRe, dt, X, Y, 1);	
  
  ## f4
  f4 = rkf(wn, w, u, v, cmatrix, invRe, dt, X, Y, 0);	
  
  ## Update
  w = 1/3*(-wn + w1 + 2*w2 + w3) + dt/6.*f4;
  w = bc(w);
	   
  wvec   = -mattovec(w);
  phivec = cmatrix\wvec;
  phi = vectomat(phivec, ngy, ngx);
  phi = bc(phi);
  
  [v u] = der2d(phi, X, Y, 'n');
  v = -v;
  u = bc(u);
  v = bc(v);
  
endfunction

function run_simulations(specs)			  

	progress_bar = waitbar(0,'Initialising ...');

	ngrid = specs.ngrid;
	lbox  = specs.lbox; 
	dx    = lbox/double(ngrid);
	dt    = specs.timestep;
	viscosity = specs.viscosity;
	nloops = specs.numloops;
	samplefreq = specs.samplefreq;
	
	x = linspace(-lbox/2, lbox/2, ngrid);
	y = x;
				
	[X, Y]    = meshgrid(x, y);
	[ngy ngx] = size(X); 
				
	cmatrix = cmat2d(X,Y);
				
	## Initial conditions
	sigma = 10;
	w1 = exp(-sigma*(X-specs.vortex1(1)).^2 - sigma*(Y-specs.vortex1(2)).^2);
	w2 = exp(-sigma*(X-specs.vortex2(1)).^2 - sigma*(Y-specs.vortex2(2)).^2);
	w3 = exp(-sigma*(X-specs.vortex3(1)).^2 - sigma*(Y-specs.vortex3(2)).^2);
	w4 = exp(-sigma*(X-specs.vortex4(1)).^2 - sigma*(Y-specs.vortex4(2)).^2);
	w  = specs.vortex1(3)*w1 + specs.vortex2(3)*w2 + specs.vortex3(3)*w3 + specs.vortex4(3)*w4;
	wm  = sum(sum(w));
	fac = 2./(dx^2)/wm;
	w   = w.*fac;

	wvec   = -mattovec(w);
	psivec = cmatrix\wvec;
	psi = vectomat(psivec, ngy, ngx);
	psi = bc(psi);
		
	[v u] = der2d(psi, X, Y, 'n'); v = -v;
	v = bc(v); u = bc(u);
	
	if specs.removefiles
			system("rm -rf *mat");
	end
	
	fplot=figure(2);
	n = 0; t=0.0;
	
	dcounter = 0;
	filename = sprintf("dump-%d.mat", dcounter);
	save(filename, '-mat-binary', "w", "psi", "u", "v", "X", "Y", "t");
	dcounter = dcounter + 1;
	
	while n <= nloops
						
		if specs.euler
			[w, u, v, psi] = eulerstep(w, u, v, cmatrix, viscosity, dt, X, Y);
		else 
			[w, u, v, psi] =  rk4step(w, u, v, cmatrix, viscosity, dt, X, Y);
		endif
				
		if rem(n, samplefreq)==0
			
			if specs.surface 
			  surf(X,Y, w)
			else
			  contour(X, Y, w, 15, 'k');
			endif
		
			filename = sprintf("dump-%d.mat", dcounter);
			save(filename, '-mat-binary', "w", "psi", "u", "v", "X", "Y", "t");
			dcounter = dcounter + 1;
		
			pause(0.001);
			waitbar(double(n)/double(nloops), progress_bar,'Integrating ... (Press Crtl-C in prompt to abort)');	
		
		endif			    
		
		t = t + dt; n++;
	endwhile
	
	waitbar(1, progress_bar,'Done');
	close(progress_bar); close(fplot);
endfunction


##
## Callback functions - note the functions must have the interface fun(guidata object, bool)
##

function get_numbgrids(obj, init=false)
	global specs;
	
	h = guidata(obj);
	
	ngrid = int64(500*get (h.grid_size, "value"));
	set (h.grid_label, "string", sprintf ("Numb. grid points: %d", ngrid));
	
	specs.ngrid = ngrid;
	
endfunction

function run_sim(obj, init=false)
	global specs
	
	run_simulations(specs);	
	
endfunction

function get_text_entry(obj, init=false)
	global specs
	
	h = guidata(obj);
	switch (gcbo)
			case {h.viscosity}
				specs.viscosity =  str2double(get(gcbo, "string")); #'value' does aply
			case {h.timestep}
				specs.timestep = str2double(get(gcbo, "string"));
			case {h.numloops}
				specs.numloops = int64(str2double(get (gcbo, "string")));
			case {h.samplefreq}
				specs.samplefreq = int64(str2double(get (gcbo, "string")));
			case {h.boxlength}
				specs.lbox = str2double(get(gcbo, "string"));
	endswitch
	
endfunction

function surface_plot(obj, init=false) 
	global specs;
	value = get(gcbo, "value");
	if value == 0
		specs.surface = false;
	elseif value == 1
		specs.surface = true;
	endif 
	
endfunction 

function remove_files(obj, init=false) 
	global specs;
	value = get(gcbo, "value");
	if value == 0
		specs.removefiles = false;
	elseif value == 1
		specs.removefiles = true;
	endif 
	
endfunction 

function set_scheme(obj, init=false) 
	global specs;
	
	h = guidata(obj);
	
	switch (gcbo)
		case {h.scheme_euler}
			specs.euler = true;
			set (h.scheme_rk4, "value", 0);
			set (h.scheme_euler, "value", 1);
		case {h.scheme_rk4}
			specs.euler = false;
			set (h.scheme_euler, "value", 0);
			set (h.scheme_rk4, "value", 1);
	endswitch
	
endfunction 

function get_vortex_entry(obj, init=false)
	global specs
	
	h = guidata(obj);
	switch (gcbo)
			case {h.init1_xcoor}
				specs.vortex1(1) = str2double(get(gcbo, "string"));	
			case {h.init1_ycoor}
				specs.vortex1(2) = str2double(get(gcbo, "string"));
			case {h.init1_strength}
				specs.vortex1(3) = str2double(get(gcbo, "string"));	
				
			case {h.init2_xcoor}
				specs.vortex2(1) = str2double(get(gcbo, "string"));	
			case {h.init2_ycoor}
				specs.vortex2(2) = str2double(get(gcbo, "string"));
			case {h.init2_strength}
				specs.vortex2(3) = str2double(get(gcbo, "string"));	
				
			case {h.init3_xcoor}
				specs.vortex3(1) = str2double(get(gcbo, "string"));	
			case {h.init3_ycoor}
				specs.vortex3(2) = str2double(get(gcbo, "string"));
			case {h.init3_strength}
				specs.vortex3(3) = str2double(get(gcbo, "string"));	
				
			case {h.init4_xcoor}
				specs.vortex4(1) = str2double(get(gcbo, "string"));	
			case {h.init4_ycoor}
				specs.vortex4(2) = str2double(get(gcbo, "string"));
			case {h.init4_strength}
				specs.vortex4(3) = str2double(get(gcbo, "string"));	
	endswitch													
				
endfunction


f=figure(1, 'position',[100,100,700,700]);

handle.grid_label = uicontrol ("style", "text",...
"units", "normalized",...
"string", "*The Merger*",...
"horizontalalignment", "center",...
"position", [0.4 0.925 0.2 0.05]);

# Main ui panel
panel_system = uipanel("title", "System setup", "position", [.1 .4 .8 .5]);
#panel_simulation = uipanel("title", "Simulation setup", "position", [.1 .4 .8 .2]);
panel_init = uipanel("title", "Initial conditions    x coord.      y coord.      strength ", ... 
"position", [.1 .15 .6 .2]);


## System setup - system panel
handle.grid_size = uicontrol ("parent", panel_system, ...
"style", "slider",...
"units", "normalized",...
"string", "Numb. grid points: ",...
"callback", @get_numbgrids,...
"value", 0.5, ...
"position", [0.05 0.75 0.4 0.05]);

handle.grid_label = uicontrol ("parent", panel_system, ...
"style", "text",...
"units", "normalized",...
"string", "Numb. grid points: 100",...
"horizontalalignment", "left",...
"position", [0.05 0.8 0.5 0.1]);

handle.boxlength_label = uicontrol ("parent", panel_system, ...
"style", "text",...
"units", "normalized",...
"string", "Box length",...
"horizontalalignment", "left",...
"position", [0.05 0.525 0.3 0.1]);

handle.boxlength = uicontrol ("parent", panel_system, ...
"style", "edit", ...
"units", "normalized", ...
"string", "Enter length", ...
"callback", @get_text_entry,...
"position", [0.2 0.53 0.2 0.075]);

handle.viscosity_label = uicontrol ("parent", panel_system, ...
"style", "text",...
"units", "normalized",...
"string", "Viscosity",...
"horizontalalignment", "left",...
"position", [0.55 0.75 0.2 0.1]);

handle.viscosity = uicontrol ("parent", panel_system, ...
"style", "edit", ...
"units", "normalized", ...
"string", "Enter viscosity", ...
"callback", @get_text_entry,...
"position", [0.7 0.75 0.2 0.075]);

handle.timestep_label = uicontrol ("parent", panel_system, ...
"style", "text",...
"units", "normalized",...
"string", "Time step",...
"horizontalalignment", "left",...
"position", [0.55 0.55 0.2 0.1]);

handle.timestep = uicontrol ("parent", panel_system, ...
"style", "edit", ...
"units", "normalized", ...
"string", "Enter time step", ...
"callback", @get_text_entry,...
"position", [0.7 0.55 0.2 0.075]);

handle.loop_label = uicontrol ("parent", panel_system, ...
"style", "text",...
"units", "normalized",...
"string", "No. loops",...
"horizontalalignment", "left",...
"position", [0.55 0.35 0.2 0.1]);

handle.numloops = uicontrol ("parent", panel_system, ...
"style", "edit", ...
"units", "normalized", ...
"string", "Enter no. loops", ...
"callback", @get_text_entry,...
"position", [0.7 0.35 0.2 0.075]);

handle.samplefreq_label = uicontrol ("parent", panel_system, ...
"style", "text",...
"units", "normalized",...
"string", "Sample freq.",...
"horizontalalignment", "left",...
"position", [0.5 0.15 0.2 0.1]);

handle.samplefreq = uicontrol ("parent", panel_system, ...
"style", "edit", ...
"units", "normalized", ...
"string", "Sample freq.", ...
"callback", @get_text_entry,...
"position", [0.7 0.15 0.2 0.075]);

handle.surfcont = uicontrol ("parent", panel_system, ...
"style", "checkbox",...
"units", "normalized",...
"string", "Surface plot",...
"value", 0,...
"callback", @surface_plot,...
"position", [0.25 0.25 0.2 0.09]);

handle.removefiles = uicontrol ("parent", panel_system, ...
"style", "checkbox",...
"units", "normalized",...
"string", "Remove files",...
"value", 0, ...
"callback", @remove_files,...
"position", [0.25 0.15 0.23 0.09]);

handle.scheme_text = uicontrol ("parent", panel_system, ...
"style", "text", ...
"units", "normalized", ...
"string", "Num. scheme", ...
"horizontalalignment", "left", ...
"position", [0.05 0.35 0.2 0.1]);

handle.scheme_euler = uicontrol ("parent", panel_system, ...
"style", "radiobutton", ...
"units", "normalized", ...
"string", "Euler", ...
"callback", @set_scheme, ...
"position", [0.05 0.25 0.1 0.1]);


handle.scheme_rk4 = uicontrol ("parent", panel_system, ...
"style", "radiobutton", ...
"units", "normalized", ...
"string", "RK4", ...
"callback", @set_scheme, ...
"position", [0.05 0.15 0.1 0.1]);

## Initial conditions
handle.init1_label = uicontrol ("parent", panel_init, ...
"style", "text",...
"units", "normalized",...
"string", "Vortex 1",...
"horizontalalignment", "left",...
"position", [0.05 0.7 0.2 0.1]);

handle.init1_xcoor = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.3 0.675 0.15 0.15]);

handle.init1_ycoor = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.5 0.675 0.15 0.15]);

handle.init1_strength = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.7 0.675 0.15 0.15]);

handle.init2_label = uicontrol ("parent", panel_init, ...
"style", "text",...
"units", "normalized",...
"string", "Vortex 2",...
"horizontalalignment", "left",...
"position", [0.05 0.5 0.2 0.1]);

handle.init2_xcoor = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.3 0.475 0.15 0.15]);

handle.init2_ycoor = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.5 0.475 0.15 0.15]);

handle.init2_strength = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.7 0.475 0.15 0.15]);

handle.init3_label = uicontrol ("parent", panel_init, ...
"style", "text",...
"units", "normalized",...
"string", "Vortex 3",...
"horizontalalignment", "left",...
"position", [0.05 0.3 0.2 0.1]);

handle.init3_xcoor = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.3 0.275 0.15 0.15]);

handle.init3_ycoor = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.5 0.275 0.15 0.15]);

handle.init3_strength = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.7 0.275 0.15 0.15]);

handle.init4_label = uicontrol ("parent", panel_init, ...
"style", "text",...
"units", "normalized",...
"string", "Vortex 4",...
"horizontalalignment", "left",...
"position", [0.05 0.1 0.2 0.1]);

handle.init4_xcoor = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.3 0.075 0.15 0.15]);

handle.init4_ycoor = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.5 0.075 0.15 0.15]);

handle.init4_strength = uicontrol ("parent", panel_init, ...
"style", "edit", ...
"units", "normalized", ...
#"string", "Sample freq.", ...
"callback", @get_vortex_entry,...
"position", [0.7 0.075 0.15 0.15]);

## Simulation button
handle.button = uicontrol ("style", "pushbutton", ...
"units", "normalized",...  # wrt the panel
"string", "Run", ...
"callback", @run_sim, ...
"position", [0.8 0.15 0.1 0.1]);

guidata (gcf, handle)

