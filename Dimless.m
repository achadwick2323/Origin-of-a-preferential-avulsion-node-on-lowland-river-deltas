%% master function    
function [Output] = Dimless(varargin)

    %% matlab settings & notes
    
        % Dimensionless backwater morphodynamics (over a single avulsion cycle)
        % This code runs by itself, and is also used as a subroutine by MultiDimless script.
        % Austin chadwick  3/18/19
        % Contact: achadwick@caltech.edu

        % A note about input arguments:
            % varargin{1} = parameters
            % varargin{2} = InitialConditions

        % Figure settings
        set(0,'DefaultFigureWindowStyle','docked')
        set(0,'DefaultLineLinewidth',1);
        set(0,'DefaultAxesFontSize',14);
        set(0,'DefaultLineMarkerSize',6);
        format long g
        warning('off')
        % colors
        red   = [0.915294117647059, 0.28156862745098, 0.287843137254902];  
        blue  = [0.34666666666666, 0.536, 0.9];
        green = [0.44156862745098, 0.749019607843137, 0.432156862745098];

        % use this command if you want to fold code
        %org.netbeans.api.editor.fold.FoldUtilities.collapseAll(org.netbeans.api.editor.fold.FoldHierarchy.get(com.mathworks. mlservices.MLEditorServices.getEditorApplication.getActiveEditor.getTextComponent()));
            
        % output filename
        folder = 'OutputFolder';
        filepath = fullfile(folder, ['output_' datestr(datetime('now'),'yymmdd-HHMMSS')]);
        
        % printing
        plotbin = 0; % plot it 
        savebin = 1; % save it

    %% set up parameters
    
        % if there are no input parameters specified, or they say "default", use these defaults
        if length(varargin) == 0 || strcmp(varargin{1},'default');
        parameters = ParameterBuilder_fun;  
        
        % if there are input parameters, they are the first input argument
        else
        parameters = varargin{1};
        end
        
        % unpack parameters structure to workspace
        sigma = parameters.sigma;
        fields = fieldnames(parameters, '-full');
        for i = 1:length(fields)
            commandLine = sprintf('%s = parameters.%s;', fields{i}, fields{i});
            eval(commandLine);
        end
        clear fields
        
        % sediment transport relation, depending on input string
            % EngelundHansen
            if strcmp(SedimentTransport.Relation,'EngelundHansen');
            SedimentTransport_fun = @(tau,hn,Cf,parameters) EngelundHansen_fun(tau,hn,Cf,parameters);
            end
            parameters.SedimentTransport.SedimentTransport_fun = SedimentTransport_fun; % write to parameters
                
        % flow width relation, depending on parameter string
            % constant
            if strcmp(Width.FlowWidthRelation,'Constant');
            FlowWidth_fun = @(x,xm,parameters) ConstantWidth_fun(x, xm, parameters); 
            % flare
            elseif strcmp(Width.FlowWidthRelation,'Flare'); 
            FlowWidth_fun = @(x,xm,parameters) FlareWidth_fun(x, xm, Width.Flare.theta_Flow, parameters);    
            end
            parameters.Width.FlowWidth_fun = FlowWidth_fun; % save function to parameters   
            
        % depo width relation, depending on parameter string
            % constant
            if strcmp(Width.DepoWidthRelation,'Constant');
            DepoWidth_fun = @(x,xm,parameters) ConstantWidth_fun(x, xm, parameters); 
            % flare
            elseif strcmp(Width.DepoWidthRelation,'Flare'); 
            DepoWidth_fun = @(x,xm,parameters) FlareWidth_fun(x, xm, Width.Flare.theta_Depo, parameters);    
            % step
            elseif strcmp(Width.DepoWidthRelation,'Step');
            DepoWidth_fun = @(x,xm,parameters) StepWidth_fun(x, xm, parameters); 
            end
            parameters.Width.DepoWidth_fun = DepoWidth_fun; % save function to parameters   

        % if flow relation is multi-phase, we need to calculate the
        % discretization of the distribution into into logarithmically spaced bins
        if strcmp(FlowRegime.Scheme,'MultiPhase');
            FlowRegime.MultiPhase = BinMulti_fun(FlowRegime.MultiPhase);
            parameters.FlowRegime.MultiPhase = FlowRegime.MultiPhase; 
        end
        
        % flow relation depending on input string, and calculating the
        % arithmetic mean einstein number
            % MultiPhase
            FlowRegime_fun = @(hn,Cf,Tprev_flow,parameters) MultiPhaseFlow_fun(hn,Cf,Tprev_flow,parameters);
                temp.hn = exp(FlowRegime.MultiPhase.psii);
                temp.qsn = SedimentTransport_fun(taun_bf*temp.hn,temp.hn,FlowRegime.MultiPhase.Cfi,parameters);
            qsn_mean = sum(temp.qsn.*FlowRegime.MultiPhase.Fi); 

            % finish up
            clear temp
            parameters.FlowRegime.FlowRegime_fun = FlowRegime_fun; 
            parameters.qsn_mean = qsn_mean;
            
            % unpack
            Sa = SedimentTransport.Sa;
            
    %% set up initial conditions
    
        % if there are no input initial conditions or they say "default", build from scratch
        if length(varargin) <= 1 || strcmp(varargin{2},'default')
            InitialConditions = Init_fun(parameters);
                
        % if there input initial conditions, use them
        else
        InitialConditions = varargin{2};
        end
               
        % update initial conditions in the parameter file 
        parameters.InitialConditions = InitialConditions; 
                                                     
    %% preparations
        
        % initialize time & time-keeping 
        t = 0;
        dt = 0;
        Tprev.flow = t; % Time since previous flow. When this quantity reaches a value of Te, then we change the discharge.
        cond_newhn = 1; % When this variable equals 1, we update the discharge. 
        % bed profile
        xs  = InitialConditions.xs; % sediment shock
        xb = InitialConditions.xb; % toe 
        etab = InitialConditions.etab; % toe elev
        eta = InitialConditions.eta; % bed elev
        etaF= InitialConditions.etaF; % floodplain elev
        xff   = InitialConditions.xff; % false floor loc
        etaff = InitialConditions.etaff; % false floor elev
        etaFff= InitialConditions.etaFff; 
            % floodplain fit
            [~,pF] = Floodplain_fun(xs*sigma,eta,1,xs,parameters); 
            
        % constrain the delta toe & false floor profile
        ind = xff > xb; % needs to be xfull, otherwise we cant add pieces to xff.  
        xff = [xb; xff(ind)]; % false floor x
        etaff = [etab; etaff(ind)]; % false floor profile
            
        % forcing
            % sea level
            xio = InitialConditions.xio; % sea level
            % normal-flow conditions
            if strcmp(FlowRegime.Scheme,'MultiPhase');
            [~,ind] = min(abs(1 - exp(FlowRegime.MultiPhase.psii)));
            hn = exp(FlowRegime.MultiPhase.psii(ind));
            Cf = FlowRegime.MultiPhase.Cfi(ind); 
            end
           
        % auxiliary things (not saved)
            % ocean shocks, where hydraulics start
            xso = xs;
            xbo = xb;
            % slopes
            S   = SlopeCalc_3pt_fun(xs*sigma,eta); 
            Sff = -diff(etaff)./diff(xff); Sff = [Sff; Sff(end)]; 
                Sb = 0;

            % widths
            xm = FindMouth_fun(xs*sigma,eta,xs,xb,etab,etaF,xio,parameters); % mouth, where flare begins          
            [w,wx] = FlowWidth_fun(xs*sigma,xm,parameters); % flow width
            ws     = DepoWidth_fun(xs*sigma,xm,parameters); % depo width
            ws_fore = nanmean(DepoWidth_fun([xs,xb],xm,parameters)); % effective depo width of the foreset (average of the width @ tf & fb breaks)
            
        % run water
        xw = xs*sigma; % grid nodes for hydraulics, including downstream of shock. 
        etaw = eta;
        h_bf = ones(size(eta));
        [h,tau,xpin,Fr] = Backwater_fun(xs*sigma,etaw,hn,xio,Cf,w,wx,xm,0,parameters); 
            hc = Backwater_fun(xs*sigma,etaw,1,xio,Cf,w,wx,xm,0,parameters); 
       
        % route sediment
        [qs,qs_feed] = SedimentTransport_fun(tau,hn,Cf,parameters);    
            qsfeedmeanobs = qs_feed;
        
        % avulsions
        etaref  = InitialConditions.etaref; % set up reference profile
            hcref   = InitialConditions.hcref; 
            etaFref = InitialConditions.etaFref; 
        criterion = nan(size(eta)); % criterion is not evaluated
            
            % initialize to-print matrices
            t_print    = nan(1,Mprint);
            xs_print   = nan(1,Mprint);
            xb_print   = nan(1,Mprint);
            xm_print   = nan(1,Mprint);
            etab_print = nan(1,Mprint);
            xib_print  = nan(1,Mprint); 
            qsfeed_print = nan(1,Mprint);
            eta_print       = nan(N,Mprint);
            h_print         = nan(N,Mprint);
            w_print = nan(N,Mprint);
            ws_print = nan(N,Mprint); 
            ws_fore_print = nan(1,Mprint);
            hc_print = nan(N,Mprint);
            etaF_print      = nan(N,Mprint);
            qs_print        = nan(N,Mprint);  
            criterion_print = nan(N,Mprint);
            casenum_print = nan(1,Mprint);
            qsfeedmeanobs_print = nan(1,Mprint); % observed cumulative mean in sediment feed 
            
            % save initial timestep
            t_print(1) = t;
            xs_print(1)    = xs;
            xb_print(1)    = xb;
            xm_print(1)    = xm;
            etab_print(1) = etab;
            xib_print(1) = xio(1);
            qsfeed_print(1) = qs_feed;
            eta_print(:,1) = eta;
            h_print(:,1)   = h;
            hc_print(:,1) = hc;
            w_print(:,1) = w(1:N);
            ws_print(:,1) = ws(1:N);
            ws_fore_print(1) = ws_fore;
            etaF_print(:,1) = etaF;
            qs_print(:,1)  = qs;   
            criterion_print(:,1) = criterion;
            qsfeedmeanobs_print(1) = qsfeedmeanobs; 
                                   
    %% step through time
    testy = 1; 
    tic
    for m = 2:Mprint; % times to print
                
        for n = 1:Mtoprint % model iterations between prints
            
            % case 1) if the discharge is new ...
            % send the shock front to the ocean & trigger shock search scheme
            if cond_newhn; cond_newhn = 0; casenum = 1;
                % recall full profile of previous timestep
                xfull = [xs*sigma;  xff];
                etafull = [eta;   etaff];
                etaFfull= [etaF; etaFff]; 
                % send the shock front to ocean (maximum x-value) & interpolate
                xs = xso; 
                xb = xbo; 
                eta   = interp1qr(xfull,etafull,xs*sigma);
                etaF = interp1qr(xfull,etaFfull,xs*sigma);  
                % if = 1, begin searching for a new shock front on the next iteration
                cond_shockdetect = 0;
                
            % case 2) if we triggered the shock-search scheme ...
            % send the shock front backward based on sediment flux criterion
            elseif cond_shockdetect; cond_shockdetect = 0; casenum = 2;
                % recall full profile of previous timestep
                xfull = [xs*sigma; xff];
                etafull = [eta; etaff];
                etaFfull= [etaF; etaFff];
                % search for a new shock front upstream of the current one
                [xs,xb] = ShockDetector_fun(xs,xb,qs,Sforw,xm,xso,xbo,eta,etab,parameters);
                eta  = interp1qr(xfull,etafull,xs*sigma);
                etaF = interp1qr(xfull,etaFfull,xs*sigma); 
                
            % case 3) if discharge is unchanged & we are not searching for a new shock
            % solve exner & advance thru time!
            else; casenum = 3;
                %if testy; keyboard; testy = 0; end
                % solve exner
                xb_prev = xb; 
                Sbo = 0; % assume toe slope is zero. 
                [xs,xb,eta,xso,xbo,dt] = SolveExnerClosed_fun(xs,xb,eta,S,Sb,qs,qs_feed,w,ws,ws_fore,xso,xbo,Sbo,parameters);
                % advance through time
                t = t + dt; 
                Tprev.flow  = Tprev.flow  + dt;   
                % update floodplain
                h_bf = ones(size(eta));
                [etaF,pF] = Floodplain_fun(xs*sigma,eta,h_bf,xs,parameters); 
                    etaFff = polyval(pF,xff); 
                
            end % end if cond
                        
            % constrain the delta toe & false floor profile
            etab  =  eta(end) - Sa*(xb - xs); 
                if xb == xbo; 
                etab = 0;
                xb = xs+(eta(N)-etab)/Sa;
                xbo = xb;
                end; 
            etaFb = etaF(end)- Sa*(xb - xs); 
            ind = xfull > xb; % note: needs to be xfull, otherwise we cant add pieces to xff.  
            xff = [xb; xfull(ind)]; % false floor x (i.e., basin floor)
            etaff = [etab; etafull(ind)]; % false floor profile
            etaFff= polyval(pF,xff); % [etaFb; etaFfull(ind)];
            
            % grid for hydraulics goes to the ocean no matter where the shock is.
            ind = (xff<=xso); % indices where we run water, but not sed
            xw = [xs*sigma; xff(ind)]; 
            etaw = [eta; etaff(ind)];  
            
            % find river mouth through interesection with sea level (except for during shock-trigger scheme)
            if casenum ~=2;
            xm = FindMouth_fun(xs*sigma,eta,xs,xb,nan,etaF,xio,parameters); % river mouth
            end
                        
            % forcing
            xio = InitialConditions.xio + dxibdt*t; % sea level
            [hn, Cf, Tprev.flow, cond_newhn] = FlowRegime_fun(hn,Cf,Tprev.flow,parameters); % normal flow conditions
                cond_newhn = 0; % with this set to zero, changin the flow regime will not trigger the case 1) scenario wit the shock moving downstream. 
                
            % calculate the fluvial & false-floor slopes
            S   = SlopeCalc_3pt_fun(xs*sigma,eta);
            Sff = -diff(etaff)./diff(xff); Sff = [Sff; nan]; % fwd difference
                Sb = (xs==xso)*0 + 1*(xs<xso); % Delta toe progrades over Slope=0 if we are in the ocean, and Slope~1 if we are upstream.
                
            % calculate widths
            parameters.temp.etaw = etaw; parameters.temp.xio = xio; % needed to calculate relative submergence
            etaFtemp = [etaF; etaFff]; 
            etaFw = etaFtemp(1:length(xw)); 
            [w,wx] = FlowWidth_fun(xw,xm,parameters); % flow width
            ws     = DepoWidth_fun(xs*sigma,xm,parameters); % depo width
            ws_fore = nanmean(DepoWidth_fun([xs,xb],xs,parameters)); % effective depo width of the foreset (average of the width @ tf & fb breaks)

            % find channel depth profile 
            hc = Backwater_fun(xw,etaw,1,xio,Cf,w,wx,xm,t,parameters); 
            
            % route water over the whole profile where x <= xso
            [h,tau,xpin,Fr] = Backwater_fun(xw,etaw,hn,xio,Cf,w,wx,xm,t,parameters);
             
            % route sediment over profile
                % info for previous timestep, needed for cumulative average calculation
                qsfeedmeanobs_prev = qsfeedmeanobs;
                t_prev = t - dt; 
                qs_feed_prev = qs_feed; 
                % route!
                [qs,qs_feed] = SedimentTransport_fun(tau(1:N),hn,Cf,parameters);
                % running average of observed sediment flux (used for mass-balance error estimate)
                if t>0; 
                qsfeedmeanobs = 1./t * (qsfeedmeanobs_prev*t_prev + 0.5*(qs_feed + qs_feed_prev)*dt); 
                end
                
            % check for delta fronts
                % did we trigger a front?
                Sforw = -diff(eta)./diff(xs*sigma); Sforw = [Sforw; nan]; % fwd diff slope on topset, used only for triggering shock detection
                [xs0] = ShockDetector_fun(xs,xb,qs,Sforw,xm,xso,xbo,eta,etab,parameters);
                if (xs-xs0) > 0; 
                cond_shockdetect = 1; 
                disp('    ->shock triggered.')
                end
                % did we abandon a front? (i.e. front and toe collide)
                hforemin = 0.1; % minimum allowable foreset wedge thickness before abandonment
                if Sa*(xb-xs)<hforemin; 
                disp('    ->abandoning shock front because front & toe merged.')
                cond_newhn = 1; 
                end
                        
            % calculate avulsion criterion
            [criterion,avulsion] = AvulsionCriterion_fun(xs,eta,InitialConditions.xs,etaref,etaFref,t,h(1:N),etaF,hc(1:N),xpin,xio,xm,parameters);
            if cond_newhn; avulsion = 0; end; % don't allow avulsions when it's the first timestep of discharge in a flow, because our grid isnt adjusted yet. 
            if avulsion; break; end; % end for-loop if we triggered an avulsion
                        
        end % end for n
        
        % break if we had an avulsion
        if avulsion; break; end
        
        % break if the river mouth is located upstream of x = 1
        if xm < 1; break; end; 
        
        % break if the basin depth drops below the basin floor
        if xio < 0; break; end;
      
        % save on occasional timesteps
        t_print(m) = t;
        xs_print(m) = xs;
        xb_print(m) = xb;
        xm_print(m) = xm;
        etab_print(m) = etab; 
        xib_print(m) = xio;
        qsfeed_print(m) = qs_feed;
        eta_print(:,m) = eta;
        h_print(:,m) = h(1:N);
        w_print(:,m) = w(1:N);
        ws_print(:,m) = ws(1:N);
        ws_fore_print(m) = ws_fore;
        hc_print(:,m) = hc(1:N);
        etaF_print(:,m) = etaF(1:N);
        qs_print(:,m) = qs;
        criterion_print(:,m) = criterion;
        casenum_print(m) = casenum;
        qsfeedmeanobs_print(m) = qsfeedmeanobs;
                
    end % end for m
    toc
            
    %% if we triggered an avulsion, find where & when
    
        % Save avulsion conditions
        if avulsion
            % find where and when
            [~,ind] = nanmax(criterion);
            xA = xs*sigma(ind);
            tA = t; 
            etaA = eta(ind); 
            % build reference functions
            etaAfun = spline_fun(xs*sigma,eta,parameters); % bed
            etaFAfun = spline_fun(xs*sigma,etaF(1:N),parameters); % levee
            hcfun = spline_fun(xs*sigma,hc(1:N),parameters); % channel depth
            % save avulsion matrix
            Avulsion.xA = xA; % avulsion location
            Avulsion.tA = tA; % avulsion time
            Avulsion.xs = xs; % shoreline coordinate
            Avulsion.xb = xb; % toe coordinate
            Avulsion.xm = xm; % mouth
            Avulsion.eta = eta; % bed elevations
            Avulsion.etaF = etaF; 
            Avulsion.etab = etab;
            Avulsion.qsfeed = qs_feed; 
                Avulsion.etaref = etaAfun;
            Avulsion.etaF = etaF; % levee elevation
                Avulsion.etaFref = etaFAfun; % function for levees at avulsion
            Avulsion.xff = xff;
            Avulsion.etaff = etaff;
            Avulsion.etaFff = etaFff;
            Avulsion.h = h(1:N); % depths
            Avulsion.w = w(1:N); 
            Avulsion.ws = ws(1:N);
            Avulsion.ws_fore = ws_fore;
            Avulsion.hc = hc(1:N);
                Avulsion.hcref = hcfun;
            Avulsion.qs = qs; % transport rates
            Avulsion.criterion = criterion; % criterion
            Avulsion.xio = xio; % water surface elev at ds end
            Avulsion.qsfeedmeanobs = qsfeedmeanobs;
            
        end % end if avulsion       
        
        % Save final conditions, regardless of whether there was an avulsion or not. 
            % build reference function s
            etafun_final = spline_fun(xs*sigma,eta,parameters); % bed
            etaFfun_final = spline_fun(xs*sigma,etaF,parameters); % levee
            hcfun_final = spline_fun(xs*sigma,hc(1:N),parameters); % channel depth
            % save final conditions matrix
            FinalConditions.t = t; % avulsion time
            FinalConditions.xs = xs; % shoreline coordinate
            FinalConditions.xb = xb; % toe coordinate
            FinalConditions.xm = xm;
            FinalConditions.eta = eta; % bed elevations
            FinalConditions.etaF = etaF; 
            FinalConditions.etab = etab;
            FinalConditions.qsfeed = qs_feed; 
                FinalConditions.etaref = etafun_final;
            FinalConditions.etaF = etaF; % levee elevation
                FinalConditions.etaFref = etaFfun_final; %  function for levees at avulsion
            FinalConditions.xff = xff;
            FinalConditions.etaff = etaff;
            FinalConditions.etaFff = etaFff;
            FinalConditions.h = h; % depths
            FinalConditions.hc = hc;
                FinalConditions.hcref = hcfun_final;
            FinalConditions.qs = qs; % transport rates
            FinalConditions.criterion = criterion; % criterion
            FinalConditions.xio = xio; % water surface elev at ds end
            FinalConditions.qsfeedmeanobs = qsfeedmeanobs;

    %% finish up

        % timing the run
        computer_runtime = toc; % [sec]

        % calculate the error in mass balance, assuming that the current profile extends to the tf break
        error.MassBalanceNoFlowVarCorrection = MassBalanceError_fun(InitialConditions,FinalConditions,1,parameters);
            
        % compile output structure
        Output.eta = eta_print;
        Output.h   = h_print;
        Output.qs  = qs_print;
        Output.xs  = xs_print;
        Output.xb  = xb_print;
        Output.xm  = xm_print;
        Output.etab= etab_print;
        Output.qsfeed = qsfeed_print; 
        Output.criterion = criterion_print;
        Output.etaF = etaF_print;
        Output.t = t_print;
        Output.xio = xib_print;
        Output.w = w_print;
        Output.ws = ws_print;
        Output.ws_fore = ws_fore_print;
        Output.hc = hc_print;
        Output.parameters = parameters;
        Output.InitialConditions = InitialConditions;
        Output.FinalConditions = FinalConditions;
        Output.casenum = casenum_print;
        Output.qsfeedmeanobs = qsfeedmeanobs_print;
        Output.error = error;
        Output.computer_runtime = computer_runtime;
        Output.filepath = filepath;
        
            % if there was an avulsion, save it
            if avulsion; Output.Avulsion = Avulsion; end
                             
        % save, if we want
        if savebin
        disp(['--> saving output to ' filepath])
        save(filepath,'-struct','Output');
        end
                                           
end
                        
%% solve backwater equation
function [h,tau,xpin,Fr] = Backwater_fun(x,eta,hn,xis,Cf,w,wx,~,t,parameters)
   
    % setup
    sigma = parameters.sigma;
    N = parameters.N;
    Nw = length(x); % number of nodes to run water thru
    taun_bf = parameters.taun_bf;
    Cf_bf= parameters.Cf_bf;
    Frn_bf = parameters.Frn_bf;
    Sthresh = parameters.SedimentTransport.Sthresh;
    
    % location of front & toe
    xs = x(N); % front
    
    % calculate slope over profile considered
    S = -diff(eta)./diff(x); S = [S; S(end)];  % forwards diff
        % upstream of foreset, use 3-pt central difference central (this is more stable when exner causes numerical diffusion)
        ind = (1:N);
        S(ind) = SlopeCalc_3pt_fun(x(ind),eta(ind)); 

    % calculate normal flow froude number for this discharge
    Frn2 = Frn_bf^2 * ( Cf_bf/Cf ); 
    
    % functions
    Fr2_fun = @(i,h) Frn2 * (w(1)./w(i)).^2 .* (hn./h).^3;  % [-] squared froude number function
    f_fun   = @(i,h) ( S(i) + Fr2_fun(i,h).*(h./w(i).*wx(i) - 1./Frn2) )./... % [-] variable-width backwater equation
                          ( 1    - Fr2_fun(i,h)     ); 
                                            
    % boundary condition @ i = N, water depth is fixed at the downstream end
    h(Nw,1) = xis - eta(Nw);
            
    % interior for i = 1:N-1, upwind-scheme backward euler
    % consider making the
    for i = fliplr(1:Nw-1);

        % backwater solution
        dx = x(i+1) - x(i);
        h_pred  = h(i+1) -     dx* f_fun(i+1,h(i+1)); % prediction for depth @i
        h(i)    = h(i+1) - 0.5*dx*(f_fun(i+1,h(i+1)) + f_fun(i,h_pred) ); % corrected prediction for depth@i
        
        % if we are at the foreset or an abandoned foreset, use borda carnot losses relation
        if (x(i) == xs); 
        iu = (i); % upstream index
        id = (i+1); % downstream index
        hu = CarnotForeset_fun(eta(id),w(id),h(id),eta(iu),w(iu),Fr2_fun(id,h(id)));
        h(i) = hu; 
        end
             
    end
       
    % calculate froude and shields number
    Fr2 = Fr2_fun((1:Nw)',h);
    Fr = sqrt(Fr2); % [-] froude number
    tau = taun_bf * h .* Fr2/Frn2; % [-] shields number
        
    % location where water surface is pinned
    xpin = x(find( h == (xis-eta), 1, 'first'));
    
    % halt simulation if h is not real. this usually means flow went
    % supercritical. save final timestep to current folder for
    % troubleshooting.
	if any(~isreal(h)); 
        FinalTimestep = ws2struct; 
        filename = [mfilename '_CrashedHydraulics_' datestr(datetime('now'),'yymmdd-HHMMSS')];
        save(filename,'-struct','FinalTimestep');
        disp('(!) Crashed in Hydraulics. Saving final timestep to:');
        disp(num2str(filename));
        keyboard;
    end
                               
end

%% solve exner equation
function [xs_new,xb_new,eta_new,xso_new,xbo_new,dt] = SolveExnerClosed_fun(xs,xb,eta,S,Sb,qs,qs_feed,w,ws,ws_fore,xso,xbo,Sbo,parameters)

    % unpack parameters
    dt_max = parameters.dt_max;
    N = parameters.N; % unpack parameters
    sigma = parameters.sigma;
    Sa = parameters.SedimentTransport.Sa;
    qsn_mean = parameters.qsn_mean; 
    s = parameters.s; % subsidence rate
    au = parameters.au;
    dsig = sigma(N) - sigma(N-1); 
    % unpack unknown variables 
    eta_new = nan(N,1);
    xs_new  = nan;
    xb_new  = nan;
    % get timestep
    dt = min([timestep_fun(xs,eta,S,qs,qs_feed,w,ws,ws_fore,xb,parameters)], dt_max); % CFL condition, or maximum
    
    % shock condition gives new t/f break
    i = N;
    gamma = -s ...
            -1./xs * 1./(ws(i)*qsn_mean) * (w(i).*qs(i) - w(i-1).*qs(i-1))./dsig;
    xsdot = 1./(Sa-S(i)) * ( 1./(ws_fore*qsn_mean) * (w(N)*qs(N))./(xb-xs) ...
                             - gamma); 
    xs_new = xs + xsdot*dt;
    
    % exner @ upstream node uses a feed rate
    i = 1; 
    eta_new(i) = eta(i) ...                      % prev elev
               - S(i)*(xs_new - xs)*sigma(i) ... % grid stretching
               - s.*dt ...                       % subsidence
               - dt./xs * 1./(ws(i)*qsn_mean) * (     au*(w(i).*qs(i) - w(i).*qs_feed  )./dsig ...
                                               + (1-au)*(w(i+1).*qs(i+1) - w(i).*qs(i))./dsig );
   
    % exner @ interior nodes use a weighted central diff scheme
    i = 2:(N-1);
    eta_new(i) = eta(i) ...                      % prev elev
               - S(i).*(xs_new - xs).*sigma(i) ... % grid stretching
               - s.*dt ...                       % subsidence
               - dt./xs * 1./(ws(i)*qsn_mean) .* (     au.*(w(i).*qs(i) - w(i-1).*qs(i-1))./dsig ...
                                                + (1-au).*(w(i+1).*qs(i+1) - w(i).*qs(i))./dsig );
    % exner @ downstream node is pure upwind in space
    i = N;
    eta_new(i) = eta(i) ...                      % prev elev
               - S(i)*(xs_new - xs)*sigma(i) ... % grid stretching
               - s.*dt ...                       % subsidence
               - dt./xs * 1./(ws(i)*qsn_mean) * (w(i).*qs(i) - w(i-1).*qs(i-1))./dsig; 

    % continuity solves for foreset-bottomset break
    dxb = 1/(Sa-Sb) * ( (eta_new(N) - eta(N)) + Sa*(xs_new - xs) );
    xb_new = xb + dxb;
        % progradation of delta toe over an abandoned foreset
        cond_overtake = (xb_new > xso) && (xb <= xso); % if we prograded across the toe.
        if cond_overtake
            dxb_max = (xso - xb); % the horizontal length of the first rectangle
            dxbo = 1/(Sa-Sbo) * ( (eta_new(N) - eta(N)) + Sa*(xs_new - xs) - (Sa-Sb)*dxb_max ); % the horizontal length of the second rectangle
            xb_new = xbo + dxbo; % new toe is beyond the previous ocean toe. 
        end
        
    % ocean shocks
    if xb_new >= xbo; % if the toe is advancing beyond the ocean toe, ocean shocks are the new shocks
        xbo_new = xb_new;
        xso_new = xs_new;
    else % if the toe is not advancing, ocean shocks unchanged.
        xbo_new = xbo;
        xso_new = xso; 
    end
    
    % if code crashes, save final timestep for troubleshooting
    if ~isreal(xs_new) || any(~isreal(eta_new)) || isnan(xs_new) || isnan(xb_new); 
        FinalTimestep = ws2struct; 
        filename = [mfilename '_CrashedExner_' datestr(datetime('now'),'yymmdd-HHMMSS')];
        save(filename,'-struct','FinalTimestep');
        disp('(!) Crashed in Exner. Saving final timestep to:');
        disp(num2str(filename));
        keyboard; 
    end; % end if
    
end

%% estimate floodplain elevation
function [etaF,pF] = Floodplain_fun(x,eta,h_bf,xs,parameters) % floodplain elevation (same rate of change as bed, but different slope in bw)

        % unpack
        N = parameters.N;
        sigma = parameters.sigma;
        
        % LLS fit to bed profile, force fit through some reference point
        x0guess = xs-1; 
        [~,ind0] = min(abs(x-x0guess));
            x0 = x(ind0);
            eta0 = eta(ind0);
        % indices to consider in the fit
        L = xs-x; % distance upstream from shock
        ind = (L < 1); % only use a certain range for the fit
        ind = ind & ~isnan(eta); 
        p = polyfitZero(x(ind)-x0,eta(ind)-eta0,1);
        etafit = eta0 + polyval(p,x-x0);
        % add depth 
        etaF = etafit + 1; % % 1; 
            % calculate eqn of this line to use for etaFff.
            pF = polyfit(x,etaF,1); 
   
        
    end

%% calculate avulsion criterion 
function [criterion,avulsion] = AvulsionCriterion_fun(xs,eta,xs_init,etaref,etaFref,t,h,etaF,hc,xpin,xio,xm,parameters)
% output
% criterion = N x 1 vector of criterion as a function of space
% avulsion = 1 if the threshold is reached, =0 if not. 

    % unpack
    sigma = parameters.sigma;
    hfill = parameters.hfill;
    N = parameters.N;
    hcref = parameters.InitialConditions.hcref; 
        
    % calculate avulsion criterion (setup)
    delta_eta = (eta+hc) - max([etaref(xs*sigma)+hcref(xs*sigma), xio*ones(size(sigma))],[],2); % superelevation of floodplain relative to to 1) reference lobe surface or 2) sea level, which ever reference is higher at each location. 
        criterion = delta_eta./hc; % normalized by local flow depth
        
    % conditions where avulsions cannot be triggered.
    indno = (xs*sigma >= xm) ... % upstream of mouth
        | (xs*sigma < 1) ... % too far upstream in model domain
        | (xs*sigma == xs*sigma(N)) ... % at topset-foreset break
        | (xs*sigma == xs*sigma(N-1)); % at node immediately upstream of topset-foreset break
    criterion(indno) = nan; % set criterion to nan if its in the "indno" zone.
    
    % do we trigger an avulsion?
    avulsion = any(criterion(~indno) >= hfill ...  % filling past the threshold
                   & t > 0.25);                    % and we are not too early in the model run (brief time to allow the step at avulsion node to erode. 
               
end

%% model parameters
function [parameters] = ParameterBuilder_fun()

    % Discretization of time & space
    N  = 200;	% [-] number of topset nodes
    sigma = linspace(0,1,N)'; % [-] dimensionless downstream distance
    nu_CFL = 0.001; % [-] CFL parameter. If this is too big, then errors in mass balance go up and instabilities may occur. 
    au = 0.75; % [-] upwind coefficient (1 for pure upwind, 0 for forward, 0.5 for central)              
    dt_max = 1e-3; % [-] max allowable timestep. Please adjust manually such that this does not exceed the flow duration parameter Te. 
    Mprint = 1e4; % [-] number of timesteps to save (assuming no avulsion triggered)
    Mtoprint = 1000; % [-] number of timesteps between prints.
    
    % Bankfull flow properties
    Cf_bf  = 0.01; % [-] friction factor for bankfull condition
    Frn_bf = 0.17; %  [-] bankfull normal-flow Froude number
    taun_bf = 1; % [-] bankfull normal flow shields number, tau*_n,bf = (hn_bf*Sn_bf) / (R*D) = hn_bf/(R*D) * Cf_bf * Frn_bf

    % Widths
    Width.FlowWidthRelation = 'Flare'; % 'Constant' or 'Flare'
    Width.DepoWidthRelation = 'Step'; % % 'Constant', 'Flare', or 'Step
        % For flare case
        Width.Flare.WidthToDepthRatio = 20; % [-] Ratio of width to depth in normal flow reac
        Width.Flare.theta_Depo = 30; % [deg] spreading angle of sediment @ mouth with flare relation % (50)
        Width.Flare.theta_Flow = 15; % 0.1; % [deg] spreading angle of water @ mouth with flare relation % (3)
        % For "Step" case, there is a step-change in depositional width at the river mouth.
        Width.Step.WidthAtMouth = 1; % [-] Depositional width downstream of mouth normalized by depositional width upstream of mouth. If equals 1, it is same as constant width case. 

    % Sediment transport
    SedimentTransport.tauc = 0.0495; % [-] shields number for threshold of motion
    SedimentTransport.Sa = 5; % [-] avalanche slope for shock capturing
    SedimentTransport.Sthresh = 5; % 0.5 * SedimentTransport.Sa; 
    SedimentTransport.Relation = 'EngelundHansen'; % 'WongParker', 'EngelundHansen'
       
    % Flow Variability
    FlowRegime.Scheme = 'MultiPhase';
        % sampling timescale
        FlowRegime.MultiPhase.Te = 1e-3; % [-] dimensionless flow duration 
        FlowRegime.MultiPhase.numbins = 20; % [-] number of discharge bins
        FlowRegime.MultiPhase.Cfi = Cf_bf*ones([FlowRegime.MultiPhase.numbins,1]); % [-] assume friction factor is the same for all discharges
        % the distribution itself
        FlowRegime.MultiPhase.F_bf = 0.05; % [-] bankfull exceedance probability
        FlowRegime.MultiPhase.psi_std = 0.5; % [-] std of log of stage height, for log-normal distribution
            FlowRegime.MultiPhase.CV = sqrt(exp(FlowRegime.MultiPhase.psi_std.^2)-1); % [-] coefficient of variation of stage height
            FlowRegime.MultiPhase.psi_mean = psi_mean_calculator_fun(FlowRegime.MultiPhase.psi_std,FlowRegime.MultiPhase.F_bf,0); % [-] mean of log of stage height, for log-normal distribution
        FlowRegime.MultiPhase.hn_max = 100; % [-] max allowable normal-flow depth 
   
    % Avulsion parameters
    hfill = 0.5; % [-] dimensionless avulsion threshold  
        upbin = 1; % Set whether we allow (upbin=0) or do not allow (upbin=1) avulsions upstream of initial shoreline. We never observed avulsions downstream of initial shoreline, so setting to 1 allows code to assume this & speed up. 
    
    % sea-level parameters
    dxibdt = 0; % Relative sea-level rise rate
    s = 0; % [-] normalized subsidence rate    

    % initial conditions
    InitialConditions.xs  = 6; % [-] initial foreset 
    InitialConditions.hs  = 1; % [-] initial river-mouth depth
    InitialConditions.xio = 2; % [-] initial sea level
    InitialConditions.etab = 0; % [-] initial elevation of delta toe
    InitialConditions.Sff = 0; % [-] slope of the false floor
    InitialConditions.L = 20; % [-] length of the whole domain
    
    % pack up
    parameters = ws2struct;

end

%% minor functions
    
    % engelund & hansen (1967) total load formulation
    function [qs,qs_feed] = EngelundHansen_fun(tau,hn,Cf,parameters)
    % note that total load = bedload + suspended load, i.e. the transport regimes that make up a significant portion of the bed. does not include wash load. 

        % unpack
        taun_bf = parameters.taun_bf;
        tauc = parameters.SedimentTransport.tauc;
    
        % calculate
        qs = (0.05./Cf) .* tau.^(5/2);
            qs( tau <= tauc ) = 0;      % no sediment is transported if shields number is below the threshold of motion

        % ghost node for sediment feed
        taun = taun_bf * hn; % shields number @ upstream end
        qs_feed = (0.05./Cf) .* taun.^(5/2); % shields number @ upstream end
            %qs_feed( taun <= tauc) = 0; 


    end
    
    % slope calculator, 3pt stencil
    function [S] = SlopeCalc_3pt_fun(x,eta)
    % this function calculates slope for a 3pt stencil
    
      % initialize
        N = length(eta); % matrix dimensions
        A = zeros(N); % stencil
        dx = x(2) - x(1); 
    
        % 3-pt stencil
        A = A + diag(-1*ones(N-1,1),-1) ...
              + diag( 1*ones(N-1,1),1); 
        bc = [-2 2; -1 0];
        A(1:2,1:2) = bc; 
        A(end-1:end, end-1:end) = -fliplr(flipud(bc));
        S = -1/2 * 1/dx * A*eta; 
    
    end
         
    % Multi-phase flow regime
    function [hn_new, Cf_new, Tprev_flow_new, cond_newhn] = MultiPhaseFlow_fun(hn,Cf,Tprev_flow,parameters)
    
        % unpack
        Te =  parameters.FlowRegime.MultiPhase.Te; % event timescale
        Fi =  parameters.FlowRegime.MultiPhase.Fi; % fraction that each flow is run
        psii= parameters.FlowRegime.MultiPhase.psii; % nat log of the flow depths in bin i
        Cfi=  parameters.FlowRegime.MultiPhase.Cfi; % friction factor for each flow depth
        hn_max = parameters.FlowRegime.MultiPhase.hn_max; % max allowable flow depth
        
        % the new flow is the same as the previous ...
        hn_new = hn; 
        Cf_new = Cf; 
        Tprev_flow_new = Tprev_flow; % time elapsed since previous flow ended
        cond_newhn = 0; % this is 1 if we have changed discharge, and zero if discharge is the same as prev time
        
        % ... unless we have reached the end of a flow, when we switch the
        % flow randomly, weighting by the fractions
        if  (Tprev_flow >= Te);
            cond_newhn = 1;
            Tprev_flow_new = Tprev_flow - Te; 
            ind = randp(Fi,1); % random bin index drawn from weights according to Fi
            hn_new = exp(psii(ind)); 
                hn_new = nanmin([hn_new, hn_max]); % don't let new value exceed the max. 
            Cf_new = Cfi(ind); 
            
        end % end if
                
    
    end
        
    % constant width calculator    
    function [w,wx] = ConstantWidth_fun(x,xm,parameters)
    % calculate width as a function of downstream distance. Can be constant
    % width, variable width, etc. Make sure you comment out all the "scenarios"
    % except for the one you want. this is faster than using switch case :)

        % scenario 1) constant-width throughout
        w  = ones(size(x)); % [m] width is constant 
        wx = zeros(size(x));  % [-] dw/dx is zero
        
    end
    
    % flaring width calculator
    function [w,wx] = FlareWidth_fun(x,xm,theta,parameters)
    
        % unpack
        phi = (parameters.Width.Flare.WidthToDepthRatio * parameters.Frn_bf.^2 * parameters.Cf_bf)^-1; % (Lb / w_n_bf) factor to relate the horizontal scales
        
        % get relative submergence
        H0 = 1; % thickness of submerged channel
        if isfield(parameters,'temp') & (length(parameters.temp.etaw)==length(x)); % only use the relative submergence if number of nodes are the same.
            eta = parameters.temp.etaw;
            xio = parameters.temp.xio;
            H = xio - eta; % no-flow depth
        else
            H = H0 + 1*(x-xm); % approximation of depth from sea-level to submerged channel bed
        end
        H(H<H0) = H0; % depth for width-averaging cannot be less than the channel depth
        dH = H-H0; % approximation of depth from sea-level to submerged floodplain. cannot be negative. 
        
        % initialize
        w  = nan(size(x));
        wx = nan(size(x)); 
                
        % for x <= xm, river width is constant
        ind = ( x <= xm );
        w(ind) = 1;
            wx(ind) = 0;
        
        % for x > xm, river width expands linearly above a depth H0 = 1,
        % and then width is weighted.
            % widths above and below submerged channel
            w0 = 1; 
            wtheta = 1 + 2*tand(theta).*(x - xm) * phi; 
            % weighted width formulation
            wtemp = (H0./H)*w0 + (dH./H).*wtheta;
            wxtemp = 2*tand(theta)*phi + H0*w0./H.^2 - (H0./H).^2*2*tand(theta)*phi; 
            % plug it in only where it applies
            ind = (x > xm);
            w(ind) = wtemp(ind);
            wx(ind) = wxtemp(ind);
            
    end
    
    % step-change width calculator
    function [w,wx] = StepWidth_fun(x,xm,parameters)
    
        % unpack
        wm = parameters.Width.Step.WidthAtMouth;      
        
        % initialize
        w  = nan(size(x));
        wx = nan(size(x)); 
                
        % for x <= xm, river width is constant
        ind = ( x <= xm );
        w(ind) = 1;
            wx(ind) = 0;
            indend = find(ind,1,'last');
            wx(indend) = Inf;
            
        % for x > x0, river width expand instantaneously
        % note: this should only be used for depositional width:
        % instanteous change in flow width will cause instabilities in b/w
        % equation.
        ind = (x > xm);
        w(ind) =  wm;
            wx(ind) = 0;
                
    end
    
    % timestep calculator
    function [dt] = timestep_fun(xs,eta,S,qs,qs_feed,w,ws,ws_fore,xb,parameters)
    
        % unpack
        Cf_bf = parameters.Cf_bf;
        taun_bf = parameters.taun_bf;
        N = parameters.N;
        sigma = parameters.sigma;
        nu_CFL = parameters.nu_CFL;
        qsn_mean = parameters.qsn_mean; 
        Sa = parameters.SedimentTransport.Sa; 
                
        % CFL condition for exner
        deta = 0.1; % [-] scale of change in bed elevation (timestep scales linearly with this)
        [dwqs,ind] = max([abs(diff(w(1:N).*qs)); w(1)*qs_feed]);
        dx = xs*(sigma(N) - sigma(N-1)); 
        dt = nu_CFL * dx * deta * ( ws(ind)*qsn_mean ./ dwqs ); 
            dt = abs(dt); 
            
        % CFL condition for shock
        dt_shock = nu_CFL * Sa * (xb-xs)^2 * (ws_fore/ws(1)) * (qsn_mean/qs(N));
            dt_shock = abs(dt_shock);
            dt = min([dt,dt_shock]); % minimum of the two CFL guesses
                    
                    
    end
    
    % relationship between mean & std for a log-normal distribution
    function [mu] = psi_mean_calculator_fun(sigma,F_bf,plotbin)

        % assuming a lognormal distribution of flow depths and set exceedence
        % probability of bankfull
        hnstar_temp = logspace(-3,0,1000);
        F_bf_fun = @(mu,sigma) 1 - trapz(hnstar_temp,pdf(makedist('lognormal','mu',mu,'sigma',sigma),hnstar_temp)); % calculates exceedence prob of bankfull for a given lognormal distribution of hnstar
       
        % solve for relatinoship between geometric mean & standard deviation of flows
            % initialize
            mu = nan(size(sigma));        
            % solver setings
            options = optimset('display','off','tolx',1e-6); 
            mu_guess = -0.5; % first guess
            % solve for every value of gstd
            for i = fliplr(1:length(sigma))
               %i % display
               resid = @(mu) sqrt((F_bf - F_bf_fun(mu,sigma(i))).^2); 
               [mu_optimum,~,exitflag(i,1)] = fsolve(resid,mu_guess,options); 
                   if exitflag(i) == 1
                   mu(i,1) = mu_optimum; % write
                   mu_guess = mu_optimum; % redo guess
                   end
            end
            % use linear interpolation to get nan values
            ind = ~isnan(mu);
            p = polyfit(sigma(ind),mu(ind),1);
            mu_fit = polyval(p,sigma);
   
            % plot if you want
            if plotbin
                fig = figure();
                for i = 1:length(sigma)
                subplot(1,2,1); hold on; pbaspect([1 1 1]);
                    plot(mu(i),sigma(i),'.','markersize',75); 
                    plot(mu_fit(i),sigma(i),'kx','markersize',10);
                    plot(mu_fit,sigma,'k-');
                        xlabel('\psi_m_e_a_n = mean(log( h_n / h_n_,_b_f ))')
                        ylabel('\psi_s_t_d = std(log( h_n / h_n_,_b_f ))')
                        legend('solution found','solution interpolated');
                        title('relationship between geometric mean & std')
                subplot(1,2,2); hold on; pbaspect([1 1 1]);
                    if ~isnan(mu(i))
                    pd = makedist('lognormal','mu',mu(i),'sigma',sigma(i));
                    pdf_temp = pdf(pd,hnstar_temp);
                    plot(hnstar_temp,pdf_temp,'-'); 
                    else
                    pd = makedist('lognormal','mu',mu_fit(i),'sigma',sigma(i));
                    pdf_temp = pdf(pd,hnstar_temp);
                    plot(hnstar_temp,pdf_temp,'--k');                     
                    end
                        %set(gca,'xscale','log')
                        xlabel('h_n / h_n_,_b_f')
                        ylabel('probability density')
                        ylim([0 5]);
                        title('distributions')
                end
            end

            % replace nans with the fit
            ind = isnan(mu);
            mu(ind) = mu_fit(ind); 


    end
        
    % find the river mouth
    function [xm] = FindMouth_fun(x,eta,xs,xb,etab,etaF,xio,parameters)
    
        % unpack
        sigma = parameters.sigma;
        N = parameters.N;
        xs_init = parameters.InitialConditions.xs;
        Cf_bf = parameters.Cf_bf;
        Cf = parameters.Cf_bf;
        Frn_bf = parameters.Frn_bf; 
        Sa = parameters.SedimentTransport.Sa;
        w2d = parameters.Width.Flare.WidthToDepthRatio;
        phi = (w2d * Frn_bf.^2 * Cf_bf)^-1; % (Lb / w_n_bf) factor to relate the horizontal scales
       
        % the intersection between floodplain & sea-level 
        xint = intersections([0 2*max(x)]',[xio xio]',x,etaF); 
        xint = max(xint); % downstream-most point of intersection
        xm = xint; 
            % if there is no intersection, put it to the shock
            if isempty(xm); xm = xs; end
     
    end
    
    % write workspace to structure
    function WStruct=ws2struct()
    
        WSVARS = evalin('caller', 'who');
        for wscon=1:size(WSVARS,1)
            thisvar=evalin('caller', WSVARS{wscon});
            THEWORKSPACE.(WSVARS{wscon})=thisvar;
        end

        WStruct=THEWORKSPACE;
        
    end

    % round to specified significant digits
    function y=roundsd(x,n,method)

if nargin < 2
	error('Not enough input arguments.')
end

if nargin > 3
	error('Too many input arguments.')
end

if ~isnumeric(x)
		error('X argument must be numeric.')
end

if ~isnumeric(n) || ~isscalar(n) || n < 0 || mod(n,1) ~= 0
	error('N argument must be a scalar positive integer.')
end

opt = {'round','floor','ceil','fix'};

if nargin < 3
	method = opt{1};
else
	if ~ischar(method) || ~ismember(opt,method)
		error('METHOD argument is invalid.')
	end
end

e = floor(log10(abs(x)) - n + 1);
og = 10.^abs(e);
y = feval(method,x./og).*og;
k = find(e<0);
if ~isempty(k)
	y(k) = feval(method,x(k).*og(k))./og(k);
end	
y(x==0) = 0;
    end

    % default initial bed builder
    function [eta,xff,etaff] = BedBuilder_fun(xs,hs,xib_init,parameters,plotbin)

        % topset
        Sn    = 1; % normal flow slope = 1 in our dimless scheme
        sigma = parameters.sigma;
        etab = parameters.InitialConditions.etab;
        etas  = xib_init - hs;   % [-] elevation of t/f break
        eta = etas + xs*Sn*(1-sigma);   % [-] bed elevations

        % foreset
        Sa = parameters.SedimentTransport.Sa;
        xb = xs + (etas - etab)./Sa;   % initial location of f/b break        
        
        % false floor, using a node spacing similar to of the t/f break
        xmax = parameters.InitialConditions.L; % thrice the initial shock location should suffice
        xff = (xb:(xs/parameters.N):xmax)'; %linspace(xb,xmax,(xmax-xb)/(xb-xs))'; 
        Sff = parameters.InitialConditions.Sff; 
        etaff = etab - Sff*(xff - xb); 
        
        if plotbin

            % define variables to plot
            x_plot   = [xs*sigma; xff]; % [m] x-positions of bed
            eta_plot = [eta; etaff];    % [m] elevations of bed
            xib_plot = linspace(xib_init,xib_init,length(x_plot))'; % [m] elevations of basin water surface
                xib_plot( xib_init < eta_plot) = nan; 

            % plot
            figure()
            h1=plot(x_plot, eta_plot, 'k-','linewidth',1); hold on
            h2=plot(x_plot, xib_plot, '-','color','b','linewidth',1);
                xlabel(' x / L_b')
                ylabel(' z / h_n');
                legend([h1 h2],'bed','water surface')
                title('initial conditions')  

        end

    end
    
    % initial conditions function
    function [InitialConditions] = Init_fun(parameters)
    
        % unpack parameters   
        sigma = parameters.sigma;
        xs_init = parameters.InitialConditions.xs;     % [m] initial x-position of shoreline
        xib_init = parameters.InitialConditions.xio;   % [-] xio* = xio / hn,bf
        hs_init = parameters.InitialConditions.hs;     % [-] initial water depth @ shoreline. (here set to high flow normal depth)
        
        % build fluvial bed profile
        [eta_init,xff_init,etaff_init] = BedBuilder_fun(xs_init,hs_init,xib_init,parameters,0); % build bed
            % reference function 
            etaref = spline_fun(xs_init*sigma,eta_init,parameters); 
            
        % t/b location
        etab_init = etaff_init(1); 
        xb_init = xff_init(1); 
            
        % find mouth & flow widths
        xm_init = xs_init; % approximation of mouth, since we do not know the floodplain elevation yet 
        [w,wx] = parameters.Width.FlowWidth_fun(xs_init*sigma,xm_init,parameters); % flow widths
        
        % run bankfull hydraulics to get floodplain
        [h_bf_init,~] = Backwater_fun(xs_init*sigma,eta_init,1,xib_init,parameters.Cf_bf,w,wx,xm_init,0,parameters); % bankfull flow depth, i.e. levee height above bed
            h_bfref = spline_fun(xs_init*sigma,h_bf_init,parameters);
        etaF_init = eta_init + 1; 
        etaFff_init = etaff_init + 1; 
            etaFref = spline_fun(xs_init*sigma,etaF_init,parameters);
            
       % initial channel depth profile.
       %hc_init = BackwaterChannel_fun(xs_init*sigma,eta_init,1,xib_init,parameters.Cf_bf,w,wx,xm_init,0,parameters);
       hc_init = h_bf_init;
            hcref = spline_fun(xs_init*sigma,hc_init,parameters);
        
       % toe
       xb_init = xff_init(1); 
       etab_init = etaff_init(1); 
            
       % save initial conditions
        InitialConditions.xs = xs_init;
        InitialConditions.xb = xb_init;
        InitialConditions.etab = etab_init;
        InitialConditions.eta = eta_init;
        InitialConditions.hc = hc_init; 
        InitialConditions.etaF= etaF_init;
        InitialConditions.etaref = etaref;
        InitialConditions.etaFref = etaFref; 
        InitialConditions.hcref = hcref; 
        InitialConditions.xff = xff_init;
        InitialConditions.etaff = etaff_init; 
        InitialConditions.etaFff= etaFff_init;
        InitialConditions.xio = xib_init;
        
    
                    
        
    end
    
    % shock detector
    function [xs_new,xb_new] = ShockDetector_fun(xs,xb,qs,S,xm,xso,xbo,eta,etab,parameters)
        
        % unpack
        sigma = parameters.sigma;
        N = parameters.N;
        Sthresh = 1*parameters.SedimentTransport.Sthresh; % 0.75 * parameters.SedimentTransport.Sa; 
        Sa = parameters.SedimentTransport.Sa;
        
        % delta x for later
        dx = xs*(sigma(2)-sigma(1)); 

        % the shock is downstream-most node where the slope is steeper than Sthresh
        % the idea is, if the slope is steeper upstream, it will trigger on the next iteration. 
        ind1 = (S >= Sthresh);
        ind2 = (S == max(S));
        indshock = find(ind1 & ind2,1,'last');
        xs_new = xs*sigma(indshock);
        xs_new = min([xs_new, xs]); % helps in case ind = empty
        
        % don't trigger shocks upstream of xs = 1.
        if xs_new < 1; 
            xs_new = xs; % go back to the old one. 
            %disp('    -> ignoring steep slope upstream of x = 1')
        end
        
        % toe is the same, except if the shock moves upstream, then toe is moved arbitrarily close
        xb_new = xb; 
        if xs_new < xs; 
            % first point moving downstream from xs_new where slope is less
            % than the threshold slope. 
            ind3 = (1:length(sigma))' > indshock; % downstream of xs_new
            ind4 = (S<Sthresh); % slope is shallow
            indb = find( ind3 & ind4, 1, 'first'); 
            xb_new = xs*sigma(indb);
            % prescribe a min value for xb
            xb_newmin = xs_new + dx; % minimum   
            xb_new = max([xb_new, xb_newmin]);
        end
        
        % however, if we are detecting a shock one node upstream of the
        % current shock, it deserves special treatment. basically, the toe
        % does not move, and the shock is determined using the (N-1) node
        % elevation relative to the toe elevation (i.e. hte new foreset
        % thickness) and the avalanche slope.
        if indshock == (N-1);
            disp('    ->triggered (N-1) shock case')
            hfore = eta(N-1)-etab;
            xb_new = xb;
            xs_new = xb - hfore./Sa;
        end
        
%         % plot if u want
%         if xs_new < xs;
%             % plot and keyboard
%             figure(); hold on; pbaspect([1 1 1]); 
%             xlim([-1 8]);
%             plot(xs*sigma,eta,'ks-');
%             plot(xs,eta(end),'kv','markerfacecolor','k');
%             plot(xb,etab,'k^','markerfacecolor','w');
%             plot(xs_new,eta(indshock),'kv','markerfacecolor','r');
%             keyboard
%         end % if plot
                                   
    end

   % borda-carnot losses relation
    function [hu] = CarnotForeset_fun(etad,wd,hd,etau,wu,Frd2)
    % see notes for 5/5/17
    
        % define
        deta = etau - etad; % foreset thickness
        ra0 = hd * (wd/wu); % ratio of downstream cross-sectional area to upstream cross-sectional area for bankfull normal upstream conditions
    
        % coeffs for quad eqn
        a = -1/hd; 
        b = 1 - (deta/hd) + Frd2; 
        c = -Frd2*ra0; 
        
        % solve quad eqn
        hu = ( -b - sqrt(b^2 - 4*a*c) )/(2*a);
        
        % check if we violated borda-carnot relation
        %if (etau + hu) > (etad + hd); keyboard; end
        
    end
    
    % non-exceedence to fractions
    function [ Fi, psii ] = Nonexceed2Fraction_fun( Fnei, psibi )

    % indices
    N = length(psibi); 
    k = (1:(N-1))';
    
    % psi values are the mean of the bounds of the psi bins
    psii(k,1) = 0.5*( psibi(k) + psibi(k+1) ); 
    
    % fractions are the differences in the fractions finer
    Fi(k,1) = abs(Fnei(k) - Fnei(k+1));

    end

    % fractions to non-exceedence
    function [ Fnei ] = Fraction2Nonexceed_fun( Fi, psibi )

        % initialize
        Fnei = zeros(size(psibi)); 

        % indices
        N = length(psibi); 

        % fractions finer
        for k = 1:(N-1)
        Fnei(1) = 1; 
        Fnei(k+1) = Fnei(k) - Fi(k);
        end

    end
        
    % measure fractional error in mass balance
    function [MassBalanceError] = MassBalanceError_fun(InitialConditions,FinalConditions,corr,parameters)
        
        % time
        t = FinalConditions.t; 
    
        % measure initial delta volume
        x1 = [InitialConditions.xs * parameters.sigma; InitialConditions.xb];
        eta1 = [InitialConditions.eta; InitialConditions.etab]; 
        V1 = trapz(x1,eta1); 
        
        % measure final delta volume
        x2 = [FinalConditions.xs * parameters.sigma; FinalConditions.xb];
        eta2 = [FinalConditions.eta; FinalConditions.etab]; 
        V2 = trapz(x2,eta2); 
        
        % measure deposit volume
        V_deposit = V2 - V1; 
        
        % measure fractional error in mass balance based on the normalized
        % difference between deposit volume and time, using the correction
        % factor "corr" that is determined by predicted vs observed flow
        % variability.
        MassBalanceError = (V_deposit - t*corr)./t*corr; 
        
    end

    % bin a multi-phase flow
    function [MultiPhase] = BinMulti_fun(MultiPhase)
    
        % if we have already prescribed psii and Fi in the parameters
        % structure, notify user and exit function
        if isfield(MultiPhase,'Fi')
           disp('    ->depths and fractions specified in parameters. distribution discretization abandoned');
           MultiPhase.psib = nan([length(MultiPhase.Fi)+1,1]); 
           MultiPhase.Fnei = nan([length(MultiPhase.Fi)+1,1]);
           return
        end
    
        % so we don't have to edit lower code that used to be in the main
        % function
        FlowRegime.MultiPhase = MultiPhase;
    
        % calculate psi bin edges and their non-exceedence fractions
        temp.pd = makedist('normal','mu',FlowRegime.MultiPhase.psi_mean,'sigma',FlowRegime.MultiPhase.psi_std);            
        temp.psi = linspace(-5,4,1000); 
        temp.Fne = cdf(temp.pd,temp.psi); % non exceedence fraction
        ind1 = find(temp.Fne < 0.01,1,'last'); 
            ind1b = find(exp(temp.psi) > 0.1,1,'first');
            ind1 = max([ind1 ind1b]);
        ind2 = find(temp.Fne > 0.995,1,'first'); 
        psib = linspace(temp.psi(ind1),temp.psi(ind2),FlowRegime.MultiPhase.numbins+1)'; % bins
        Fnei = cdf(temp.pd,psib); % non-exceedence fractions for each bin
        
        % for each bin i, calculate the discharge class magnitude (psii) and fraction (Fi) in each, and the friction factor Cfi
        [Fi, psii] = Nonexceed2Fraction_fun(Fnei,psib); % fractions and psi values for each bin
            Fi = Fi./sum(Fi); % correct Fi to make sure it sums to 1
            %Cfi = Cf_bf*ones(size(psii)); % friction factor for each flow depth
        
        % write to flow regime structure
        FlowRegime.MultiPhase.psib = psib;
        FlowRegime.MultiPhase.Fnei = Fnei;
        FlowRegime.MultiPhase.psii = psii;
        FlowRegime.MultiPhase.Fi = Fi;
        
        % write to output
        MultiPhase = FlowRegime.MultiPhase;
                
%                 % plot discretization
%                 plot(temp.psi,temp.Fne,'ks-'); hold on
%                 for i = 1:length(psib);
%                 plot(psib(i)*[1 1],ylim,'k--','linewidth',0.5);
%                 end
%                 for i = 1:length(psii)
%                 plot(psii(i)*[1 1],ylim,'k-','linewidth',2)
%                 end
%                         xlabel('\psi = ln(H_n/H_0) [-]')
%                         ylabel('fraction non-exceedence, f_n_e [-]')
%                         keyboard   
    
    end
    
    % reference spline, with or without avulsions ds of initial shoreline
    function [yref] = spline_fun(x,y,parameters)
    
        % spline
        pp = spline(x,y); 
        
        % if we do NOT allow avulsions upstream of initial shoreline
        if parameters.upbin == 1;
        pp.coefs(end,:) = nan; % do not extrapolate downstream of the spline data(note this prevents avulsoins on the final node, as well)
        yref = @(x) ppval(pp,x); % reference bed function for avulsion criterion needs to allow for interpolation beyond xs_init if youw ant to pass on in multi dimless
    
        elseif parameters.upbin == 0;
        %xs_init = parameters.InitialConditions.xs;
        %yref = @(x) (x <= xs_init).*ppval(pp,x) ... % interpolated bed
        %          + (x > xs_init) .*ppval(pp(end),xs_init); % beyond initial shoreline, reference bed = initial shoreline bed
        yref = @(x) ppval(pp,x); % reference bed function for avulsion criterion needs to allow for interpolation beyond xs_init if youw ant to pass on in multi dimless
        
        end
        
        
    end
    