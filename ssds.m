classdef ssds
    % Spectral Sensing DataSet class
    % modified from BasicClass
    % Chang Hsiung, June 2, 2017
    properties
        %Value
        nsamp
        nsampT
        nsampP
        nvar
        ncls
        nclsT
        nclsP
        LAT
        pathfname_AT
        type_Model
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = ssds(varargin) % By adding this constructor to the class definition, you can create an object in one step:
            
            if nargin == 1 
                if ischar(varargin{1})
                LAT=load(varargin{1});
                obj.pathfname_AT=varargin{1};
                elseif isa(varargin{1},'struct')
                LAT= varargin{1};  
                obj.pathfname_AT='';
                else
                    error('datatype for input not supported')
                end
                obj.nsamp=length(LAT.Atrainpk(:,1));
                obj.nsampT=length(LAT.Atrainpk(:,1));
                
                try
                obj.nsampP=length(LAT.Apred(:,1));
                catch
                obj.nsampP=NaN;    
                end
                %%%%%%
                try
                    try
                        saConc_T=LAT.PLS.Tset.saConc;
                    catch
                        saConc_T=LAT.saConc;
                    end
                catch
                    saConc_T='';
                end
                
                try
                    saConc_P=LAT.PLS.Pset.saConc;
                catch
                    saConc_P='';
                end
                %%%checking size of T and P and assign 
                if length(saConc_T)>0 & length(saConc_T)==obj.nsampT
                    obj.type_Model='PLS';
                    % check whether Atrainpk matched with that inside saConc
                    
                    %                 LAT.Atrainpk(1,:)
                    %                 LAT.saConc(1).Atrainpk(1,:)
                    
                    try
                        Atrainpk_from_saConc= cat(1,LAT.saConc.Atrainpk);
                        disp('check Atrainpk')
                        if ~isSAME_2Matrix(Atrainpk_from_saConc,LAT.Atrainpk)
                            Speak_mk('mismatch between "A" train peak vs that in structure array');
                            Speak_mk('"A" train peak inside structure array will be replaced by external "A" train peak');
%                             [LAT.saConc.Atrainpk]=mat2cell_CH_4SAinsert(LAT.Atrainpk,'row');
                            LAT.saConc=Atrainpk2saConc(LAT.Atrainpk,LAT.saConc);
                            
                            
                        end
                    catch
                        Atrainpk_from_saConc= cat(1,LAT.PLS.Tset.saConc.Atrainpk);
                        disp('check Atrainpk')
                        if ~isSAME_2Matrix(Atrainpk_from_saConc,LAT.Atrainpk)
                            Speak_mk('mismatch between "A" train peak vs that in structure array')
                            Speak_mk('"A" train peak inside structure array will be replaced by external "A" train peak');
                            [LAT.PLS.Tset.saConc.Atrainpk]=mat2cell_CH_4SAinsert(LAT.Atrainpk,'row');
                        end
                    end
                    
                    try
                        Apred_from_saConc= cat(1,LAT.PLS.Pset.saConc.Atrainpk);
                         disp('check Apred')
                        if ~isSAME_2Matrix(Apred_from_saConc,LAT.Apred)
                            Speak_mk('mismatch between "A" pred vs that in structure array');
                             Speak_mk('"A" train peak inside structure array will be replaced by external "A" pred');
                            [LAT.PLS.Pset.saConc.Atrainpk]=mat2cell_CH_4SAinsert(LAT.Apred,'row');
                        end
                    
                    
                    end
                    
                    
                    
                    
                else
                    obj.type_Model='Clsfr';
                end
                 if length(saConc_P)>0 & length(saConc_P)~=obj.nsampP
                     error('mismatch between size of PLS.Pset.saConc vs nsampP')
                 end
                %%%%%%
                obj.nvar=length(LAT.Atrainpk(1,:));
                % deal with Label Combination data structure in multi-label SVM
                if isfield(LAT,'clistclslabel_LC') && isfield(LAT,'AclassinfoP')
                    % deal with Label Combination data structure in multi-label SVM
                  obj.ncls=length(LAT.clistclslabel);  
                  obj.nclsT=length(LAT.clistclslabel);  
                  obj.nclsP=length(LAT.clistclslabel_LC); 
                else
                obj.ncls=length(LAT.clistclslabel);
                obj.nclsT=length(LAT.clistclslabel);  
                obj.nclsP=NaN; 
                end
                %%%%%%%%
                obj.LAT=LAT;
                
                
                
            else
                error('not ready to handle this case yet')
                
            end
        end
        %%%%%%%%%%%%%%%%
        
        function r = roundOff(obj)
            r = round([obj.Value],2);
        end
        
        %%%%%%%%%%%%%%%%
        function r = multiplyBy(obj,n)
            r = [obj.Value] * n;
        end
        %%%%%%%%%%%%%%%%
        %Here is an overload of the MATLAB plus function. It append Tset of obj2 to end of obj1
        function r = plus(o1,o2)  %Here is an overload of the MATLAB plus function. It defines addition for this class as adding the property values:
            %             r = o1.nsamp + o2.nsamp;
            
            disp('append obj2 to obj1')
            % check whether o1 and o2 match in "type_Model"
            if strcmp(o1.type_Model,o2.type_Model)
                type_Model=o1.type_Model;
                try
                    try
                        saConc_T_1=o1.LAT.PLS.Tset.saConc;
                    catch
                        saConc_T_1=o1.LAT.saConc;
                    end
                catch
                    saConc_T_1='';
                end
                
                try
                    try
                        saConc_T_2=o2.LAT.PLS.Tset.saConc;
                    catch
                        saConc_T_2=o2.LAT.saConc;
                    end
                catch
                    saConc_T_2='';
                end
            else
                error('obj1 and obj2 have different "type_Model"');
                
            end
%             if ~isempty(saConc_T_1) && ~isempty(saConc_T_2)
%                 type_Model='PLS';
%             elseif isempty(saConc_T_1) && isempty(saConc_T_2)
%                 type_Model='Clsfr';
%             else
%                 error('can not determine whether it is PLS or Clsfr')
%             end
            
            
            %check if nvar and clistclslabel match
            if o1.nvar==o2.nvar && ( strcmp(type_Model,'PLS') ||  ( strcmp(type_Model,'Clsfr') && isSAME_two_cstr(o1.LAT.clistclslabel,o2.LAT.clistclslabel)  )  )
                LAT_new=o1.LAT;
                LAT_new.Atrainpk=[o1.LAT.Atrainpk;o2.LAT.Atrainpk];
                LAT_new.AclassinfoT=[o1.LAT.AclassinfoT;o2.LAT.AclassinfoT];
                LAT_new.AclabelT=[o1.LAT.AclabelT;o2.LAT.AclabelT];
                %%%%%%%%%%%%%%%%%%%%%%%
                % appending saConc (only works on Tset of o2 to Tset of o1)
                if strcmp(type_Model,'PLS')
                    if isfield(LAT_new,'PLS')
                        LAT_new.PLS.Tset.saConc=[LAT_new.PLS.Tset.saConc;saConc_T_2];
                    else
                        LAT_new.saConc=[LAT_new.saConc;saConc_T_2];
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%
                try
                    RawSpectra1= o1.LAT.RawSpectra.Tset;
                catch
                    try
                        RawSpectra1= o1.LAT.RawSpectra;
                    catch
                        RawSpectra1='';
                    end
                    
                end
                
                try
                    RawSpectra2= o2.LAT.RawSpectra.Tset;
                catch
                    try
                        RawSpectra2= o2.LAT.RawSpectra;
                    catch
                        RawSpectra2='';
                    end
                    
                end
                
                if ~isempty(RawSpectra1) && ~isempty(RawSpectra2)
                    RawSpectra_new=[RawSpectra1;RawSpectra2];
                end
                
                try
                    if isa(LAT_new.RawSpectra,'struct')
                        LAT_new.RawSpectra.Tset=RawSpectra_new;
                    else
                        LAT_new.RawSpectra=RawSpectra_new;
                    end
                end
                
                obj_new=ssds(LAT_new);
                r=obj_new;
                
                
            else
                warning('mismatch in nvar or clistclslabel (only for Clsfr)');
                disp('output is still obj1')
                r=o1;
                
            end
            
            
        end  % end of "plus" method
        %%%%%%%%%%%%%%%%
        function out_obj=saveAT(obj,inp)
            %         disp('calculate nsamp, ncls, and nvar, then save as Atrainpketc_~.mat');
            %%%%%%%
            try
                corename=['_',inp.corename];
            catch
                corename='';
            end
            %%%%%%%
            snsamp=['_nsamp',num2str(obj.nsamp)];
            if ~isnan(obj.nsampP)
                snsampP=['_nsampP',num2str(obj.nsampP)];
                snsampT=['_nsampT',num2str(obj.nsamp)];
                snsamp='';
                
            else
                snsampP='';
                snsampT='';
            end
            
            if strcmp(obj.type_Model,'Clsfr')
                % deal with Label Combination data structure in multi-label SVM
                if ~isnan(obj.nclsP)
                 % deal with Label Combination data structure in multi-label SVM   
                 sncls=['_nclsT',num2str(obj.ncls),'_nclsP',num2str(obj.nclsP)];   % deal with Label Combination data structure in multi-label SVM
                else
                sncls=['_ncls',num2str(obj.ncls)];
                end
                
                
            else
                sncls='';
            end
            
            svar=['_nvar',num2str(obj.nvar)];
            %%%%
            if isempty(corename)&& ~isempty(inp.pathfname_AT)
                corename_tmp=fileparts_name_ext(inp.pathfname_AT);
                corename_tmp=find_keyword_between_markers(corename_tmp,'Atrainpketc_','.mat');
                  fname_AT=['Atrainpketc_',corename_tmp,'.mat']; 
            else
                fname_AT=['Atrainpketc_',corename,svar,sncls,snsamp,snsampT,snsampP,'.mat'];
                fname_AT=strrep(fname_AT,'__','_');
            end
            %%%%
            SAT=obj.LAT;
            save(fname_AT,'-struct','SAT');
            disp([pwd,'\',fname_AT,' has been saved !']);
            obj.pathfname_AT=[pwd,'\',fname_AT];  % update property --> pathfname_AT
            out_obj=obj;
        end %end of "saveAT" method
        %%%%%%%%%%%%%%%%%%
        % this is equivalent to the method of "gt" or ">" below
        function oTP=TP_pair(o1,o2)
            % form TP pair by using first obj as Tset and 2nd as Pset
            % see for example: test_ssds_TP_pair_method()
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find skeyP from fnLP in Atrainpk_merge_Apred()
            % in current setting, try two possibilities: 
            % "(T-" vs ")"   and   "T-" vs "_"
            inp.cmk1={'(T-','T-'};% this for extracting Pset info from "o2"
            inp.cmk2={')',  '_'};% this for extracting Pset info from "o2"
            % find skeyP from fnLP in Atrainpk_merge_Apred()
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if ~isempty(o1.pathfname_AT)&&~isempty(o2.pathfname_AT)
%               fnLT= o1.pathfname_AT;
%               fnLP= o2.pathfname_AT;
%             else
            fnLT=o1.LAT;  % always use struct (because sometimes o1.pathfname_AT may already been deleted)
            fnLP=o2.LAT;  % always use struct (because sometimes o2.pathfname_AT may already been deleted)
%             end
            
            [fnLTLP SfnLTLP]=Atrainpk_merge_Apred(fnLT,fnLP,inp);
            pathfnameLTLP=which(fnLTLP);%use this such that pathfname_AT can be filled with actual value
            oTP=ssds(pathfnameLTLP);
        end % end of "TP_pair" method
        %%%%%%%%%%%%%%%%%%
        %Here is an overload of the MATLAB gt (or ">" )function. It pairs Tset of obj2 as Pset to obj1
        % this is equivalent to the method of TP_pair
        function oTP=gt(o1,o2)
            % form TP pair by using first obj as Tset and 2nd as Pset
            % see for example: test_ssds_TP_pair_method()
            oTP=TP_pair(o1,o2);
        end % end of "gt" or ">" method
        
         %%%%%%%%%%%%%%%%%%
        function oTP=TP_pair_FalsePos(o1,o2)
            % form TP pair by using first obj as Tset and 2nd as Pset
            % where 2nd AT or o2 will be converted to all NaN cls_P
            % see for example: test_ssds_TP_pair_FalsePos_method()
            % see also: '2AT->TP_FalsePos' in ATop
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % very important: o1 and o2 should be ssds obj
            % see for example: test_ssds_TP_pair_FalsePos_method()
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isa(o1,'ssds') & isa(o2,'ssds')
                try
                    fnLT=o1.pathfname_AT;
                    fnLP=o2.pathfname_AT;
                    sT=ssds(fnLT);
                    sP=ssds(fnLP);
                catch
                    sT=o1;sP=o2;
                end
            else
                error('inputs to TP_pair_FalsePos MUST be ssds obj')
            end
                
                sP.LAT.AclassinfoT=repmat(NaN,size(sP.LAT.AclassinfoT));
                sTP=sT>sP;
                delete(sTP.pathfname_AT);
                try
                strT=find_keyword_between_markers_wlistRHS( find_keyword_between_markers(fileparts_name_ext(fnLT),'{','}'),'T-',{'_',''});;
                catch
                strT='';    
                end
                try
                strP=find_keyword_between_markers_wlistRHS( find_keyword_between_markers(fileparts_name_ext(fnLP),'{','}'),'T-',{'_',''});;
                catch
                strP='';    
                end
                inp.corename=['{','T-',strT,'_','P-',strP,'}'];
                new_sd =sTP.saveAT(inp);
                
                pathfname_AT_new= strrep(new_sd.pathfname_AT,'.mat','_clsP-NaN.mat');
                copyfile( new_sd.pathfname_AT,pathfname_AT_new);
                disp(['only this --> ',pathfname_AT_new,' has been created'])
                delete(new_sd.pathfname_AT);
            %%%%%%%%%%%%%%%
            oTP=ssds(pathfname_AT_new);
        end % end of "TP_pair_FalsePos" method
        %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%
        % remove Tset and replace it by Pset
        function out_obj=P2T(obj)
            % remove Tset and replace it by Pset
            % see for example: ...
           LATnew.Atrainpk=obj.LAT.Apred;
           LATnew.AclassinfoT=obj.LAT.AclassinfoP;
           LATnew.AclabelT=obj.LAT.AclabelP;
            LATnew.clistclslabel=obj.LAT.clistclslabel;
            try
             LATnew.wvl_standardize=obj.LAT.wvl_standardize;   
            end
            try
             LATnew.RawSpectra=obj.LAT.RawSpectra.Pset;   
            end
            try
             LATnew.saConc=obj.LAT.PLS.Pset.saConc; 
%             catch
%                 error('saConc not available for Pset')
            end
            out_obj=ssds(LATnew);
        end % end of P2T method
        %%%%%%%%%%%%%%%%%%
        % replace Atrainpk into saConc and replace external Atrainpk too
        function out_obj=Atrainpk_replace(obj,Atrainpk_new)
            % replace Atrainpk into saConc and replace external Atrainpk too
             % % do not forget to reassign obj to LHS of "=" when calling this method from outside !!!
            % see for example: Atrainpk2saConc()
            LATnew=obj.LAT;
            LATnew.Atrainpk=Atrainpk_new;
            try
                LATnew.saConc=Atrainpk2saConc(Atrainpk_new,LATnew.saConc);
            catch
                LATnew.PLS.Tset.saConc=Atrainpk2saConc(Atrainpk_new,LATnew.PLS.Tset.saConc);
            end
            out_obj=ssds(LATnew);
        end % end of Atrainpk_replace method
        %%%%%%%%%%%%%%%%%%
        % replace RawSpectra into AT objects
        function out_obj=RawSpectra_replace(obj,RawSpectra_new)
            % replace RawSpectra into AT objects
            % % do not forget to reassign obj to LHS of "=" when calling this method from outside !!!
            
            LATnew=obj.LAT;
            try
            LATnew.RawSpectra.Tset=RawSpectra_new;
            catch
            LATnew.RawSpectra=RawSpectra_new;   
            end
            out_obj=ssds(LATnew);
        end % end of Atrainpk_replace method
        %%%%%%%%%%%%%%%%%%
        % show RawSpectra and Atrainpk and allow interactively pickCurve etc
        function diagnose_AT(obj,inp)
            
            if exist('inp','var')
                
                if ischar(inp) && ( strcmp(lower(inp),'atrainpk')| strcmp(lower(inp),'rawspectra'))
                    inptmp=inp;clear inp;
                    inp.Spectra_Type=inptmp;
                    
                    % set to activate --> activate_PickSpectra_GUI_yes or activate_Clsname_CmpSpectra_gui_yes
                    setup_ShowLabel_findclosestCurve_RS_AT_gui(obj.pathfname_AT,inp);
                elseif isa(inp,'struct') && isfield(inp,'Spectra_Type') && ~isempty(inp.Spectra_Type)
                    
                    % set to activate --> activate_PickSpectra_GUI_yes or activate_Clsname_CmpSpectra_gui_yes
                    setup_ShowLabel_findclosestCurve_RS_AT_gui(obj.pathfname_AT,inp);
                    
                else
                    error('data type for "inp" not supported')
                end
                %             pathfname_AT='C:\work\JDSU\CUSTOMERS\GreenWall\ATetc\Atrainpketc_GreenWall_nvar121_ncls6_nsamp252_pp1-1stDerSGw5_pp2-SNV.mat';
                
            else
                %default
                inp.Spectra_Type='Atrainpk';
                % set to activate --> activate_PickSpectra_GUI_yes or activate_Clsname_CmpSpectra_gui_yes
               setup_ShowLabel_findclosestCurve_RS_AT_gui(obj.pathfname_AT,inp); 
                
                
            end


        end %% end of diagnose_AT() method
        %%%%%%%%%%%%%%%%%%
        function diagnose_Clsname_PushButton(obj,inp)
            % To set color and width of inner line for Pset when picked see below:
            % function Clsname_CmpSpectra_gui_4ssds() line --> set(hp_spectra_pushed_P_alt,'linewidth',0.8,'color',[0.8 0.8 0.8]);
            %
            % see also setup_ShowLabel_findclosestCurve_RS_AT_gui
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if false
                % with Apred, test with ssds
                inp.Spectra_Type='Atrainpk';
                pathfname_AT='C:\work\JDSU\CUSTOMERS\GreenWall\AT_Odd-Even_bef_rm2OLs\Atrainpketc_wRawSpectra_T-odd_P-even_GreenWall_nvar121_ncls6_pp1-1stDerSGw5_pp2-SNV_nsampT126__nsampP126_TP.mat';
                sd=ssds(pathfname_AT);
                sd.diagnose_Clsname_PushButton(inp);
                %%%%%%%%%%%%%%%%%%%%
                %inp.Spectra_Type='Atrainpk';
                inp.Spectra_Type='RawSpectra';
                pathfname_AT='C:\work\JDSU\CUSTOMERS_OSP\BMS_5P\Atrainpketc__T1130_P1201_nvar58_ncls5__pp1-1stDerSGFL7[PO2]_pp2-SNV_nsampT51_nsampP49.mat';
                sd=ssds(pathfname_AT);
                sd.diagnose_Clsname_PushButton(inp);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            inp.pick_method='Clsname_PushButton';% this setting should stay the same
            setup_ShowLabel_findclosestCurve_RS_AT_gui(obj.pathfname_AT,inp);
            
        end  %end of diagnose_Clsname_PushButton
        %%%%%%%%%%%%%%%%%%
        % apply preprocessing or pretreatment
        % for example
        % inp.PP_methods.pp1='1stDerSGw15' ;  %   'SG' 'SGw5'   '1stDerSGw17' '1stDerSGDiederick' '1stDerSGw5' 'none' '1stDer'   '2ndDer'  'SNV'
        % inp.PP_methods.pp2='SNV';   % 'SampMncn'  'none' '1stDer'   '2ndDer'  'SNV'
        % inp.corename=corename;
        % sd0=ssds(LAT);
        % out_obj=apply_PP(sd0,inp);
        
        function out_obj=apply_PP(obj,inp)
            try
                inp.pathfname_AT=['Atrainpketc_',inp.corename,'.mat'];
            catch
                inp.pathfname_AT=['Atrainpketc_','some_corename','.mat'];
            end
            outPP=pretreat_preprocess_RawSpectra_pp1_pp2_AT_TP(obj.LAT,inp);
            
            LAT_PPd=load(outPP.fname_new);
            delete(outPP.fname_new);
            inp.corename=find_keyword_between_markers(outPP.fname_new,'Atrainpketc_','.mat');
            
            sd1=ssds(LAT_PPd);
            out_obj_saveAT= saveAT(sd1,inp);
            out_obj=sd1;
            out_obj.pathfname_AT=out_obj_saveAT.pathfname_AT;
            
        end % end of apply_PP() method
        %%%%%%%%%%%%%%%%%%
        % parse Atrainpk based on sub-class info in AclabelT
         function [out]=parse_AclabelT_subcls(obj,inp)
        % after parsing, all Atrainpk samples used exactly same number of times as Pset in all TP pairs
        % need to provide the following:
        %  inp.smk1='_P';inp.smk2='';inp.Nsubcls4T=2;
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % very important --> every class Must have same number of "_P",
        %  see prep_Muscle.m
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % modified from Atrainpk_Split_Odd_Even_physical_diff_samp(pathfname_AT)
        % see also Atrainpk_parse_AclabelT_subcls(), test_ssds_parse_AclabelT_subcls_method(), ATop.m
        % see also prep_ResinKits_Molecules_J_rk_4_5_6_3N1 parse_AclabelT_subcls_PDS prep_Muscle and parse_physically_different_samples
        % see doc_Hsiung_jdsu() and prep_Muscle()
        disp('work on parsing by AclabelT_subcls for classification apps')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            switch inp.parse_method
                case 'OnePDS_EachCls'
                    out=Atrainpk_parse_AclabelT_subcls(obj.LAT,inp); % see also prep_ResinKits_Molecules_J_rk_4_5_6_3N1
                case 'OnePDS_AllCls'
                    out=Atrainpk_parse_AclabelT_subcls_OnePDS_AllCls(obj.LAT,inp);
            end
        catch
            %%% for original non-PDS parsing applications, e.g. ResinKits, see for example: prep_RK_1_3_VS1120_50Polymers_wRK5_RK6
            % prep_RK_1_3_VS1120_50Polymers_wRK5_RK6() generate 3 sets of parsed results for rk1-4, rk5-6, and rk1-6
            out=Atrainpk_parse_AclabelT_subcls(obj.LAT,inp);% for original non-PDS parsing applications, e.g. ResinKits, see for example: prep_RK_1_3_VS1120_50Polymers_wRK5_RK6
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % where out.oaTP --> object array for TP pairs
        %  and  out.tmpfolder4Save  --> folder to store all TP pairs
        
        end % end of parse_AclabelT_subcls()
        %%%%%%%%%%%%%%%%%%
        % parse Atrainpk based on sub-class info in AclabelT and leave one PDS (Physically Different Sample) out each time from each class
        % PDS: Physically Different Sample
         function [out]=parse_AclabelT_subcls_OnePDS_EachCls(obj,inp)
        % will parse PDS as one PDS in Pset for each class : inp.parse_method='OnePDS_EachCls';
        % main fuction called by this method is : parse_physically_different_samples()     
        % after parsing, all Atrainpk samples used exactly same number of times as Pset in all TP pairs
        % need to provide the following:
        %   inp.smk1='_PDS'; % this can be any length, e.g."_P" or "_PDS"
        %   inp.smk2=''; % this should not be changed to others
        %
        % try to follow these procedures in the following (see prep_Muscle):
        % sd=ssds(LAT);
        %
        % inp.corename=[fileparts_name_wo_ext(pathfname_rawdata)];;
        % sd_PPd=apply_PP(sd,inp);
        %
        %  inp.smk1='_PDS';inp.smk2=''; % inp.smk1 will be used as smk4PDS in parse_physically_different_samples
        %  sd_PPd.parse_AclabelT_subcls_PDS(inp);
        % see also prep_ResinKits_Molecules_J_rk_4_5_6_3N1(pathfname)
        % see also prep_Muscle and parse_physically_different_samples
        % see also parse_AclabelT_subcls, Atrainpk_parse_AclabelT_subcls(), test_ssds_parse_AclabelT_subcls_method(), ATop.m
        % see doc_Hsiung_jdsu() and prep_Muscle()
        disp('work on parsing by AclabelT_subcls based on leave one PDS out as Pset each time')
        inp.parse_method='OnePDS_EachCls';
        
        if ~exist('inp','var') || ~isfield(inp,'parse_method_sub')
        inp.parse_method_sub='AllPerm'; % or 'SameFQallCls' and 'AllPerm' --> default setting that will generate Nfold^ncls pairs of TP files
        end
         
        out= parse_physically_different_samples(obj.pathfname_AT,inp);% see also prep_ResinKits_Molecules_J_rk_4_5_6_3N1(pathfname)
       % out=Atrainpk_parse_AclabelT_subcls(obj.LAT,inp);
        % where out.oaTP --> object array for TP pairs
        %  and  out.tmpfolder4Save  --> folder to store all TP pairs
        end % end of parse_AclabelT_subcls_PDS()
        %%%%%%%%%%%%%%%%%%
        % parse Atrainpk based on sub-class info in AclabelT
         function [out]=parse_AclabelT_subcls_OnePDS_AllCls(obj,inp)
        % will parse PDS as one PDS in Pset for All classes : inp.parse_method='OnePDS_AllCls';     
        % after parsing, all Atrainpk samples used exactly same number of times as Pset in all TP pairs
        % need to provide the following:
        %  inp.smk1='_PDS';inp.smk2='';
        % modified from Atrainpk_Split_Odd_Even_physical_diff_samp(pathfname_AT)
        % see also Atrainpk_parse_AclabelT_subcls(), test_ssds_parse_AclabelT_subcls_method(), ATop.m
        % see also parse_AclabelT_subcls_PDS prep_Muscle and parse_physically_different_samples
        % see doc_Hsiung_jdsu() and prep_Muscle()
        disp('based on "_OnePDS_AllCls" to work on parsing by AclabelT_subcls for classification apps');
        inp.parse_method='OnePDS_AllCls';
        out=parse_physically_different_samples(obj.pathfname_AT,inp);
        % where out.oaTP --> object array for TP pairs
        %  and  out.tmpfolder4Save  --> folder to store all TP pairs
        
        end % end of parse_AclabelT_subcls()
        %%%%%%%%%%%%%%%%%%
        % merge or remove or extract certain classes in Atrainpk etc 
         function [out_obj]=merge_rm_extract_class(obj,inp)
        % after merge/removal/extract, all AclassinfoT should not have any skipping in cls seq number
        % need to provide the following:
        %  inp.cls_pick={};
        %  inp.action: 
        %                    'merge'
        %                         merged_ClsName--> 'all' --> %default
        %                         merged_ClsName--> '1st_occur' --> use the first class's name as new merged class's name
        %                    '' or 'rm' or 'remove' --> remove these classes 
        %                    'extract' --> extract these classes
        % modified from Atrainpk_merge_classes_ATop()
        % see also Run_ssds_merge_rm_extract_class --> example for running this method
        % see also Atrainpk_merge_classes_ATop asmc_Global_extract_Local,  Atrainpk_merge_classes_ATop(), Atrainpk_merge_classes_Nclusters_ATop() , ATop.m
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for running two actions see for example: prep_RK_1_3_VS1120_50Polymers(pathfname_VS,pathfname_prev_AT)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('work on merge/removal/extract certain classes in Atrainpk etc for classification apps')
        clistcls_tobe_merged=inp.cls_pick;
        out=Atrainpk_merge_classes_ATop(obj.pathfname_AT,clistcls_tobe_merged,inp);
        sd1=ssds(out.LAT);
        out_obj=sd1;
        out_obj.pathfname_AT=out.pathfname;
        
        end % end of merge_rm_extract_class()
        %%%%%%%%%%%%%%%%%%
        % extract or trim down to subset of Pset in AT files with _TP etc 
         function [out_obj]=extract_Pset(obj,inp)
             % inp.loc --> locations of Pset to be extracted
             % see also test_ssds_extract_Pset_method strrep_keyword_between_markers_wlistRHS
             disp('extract certain samples out of Pset');
             locTrim=setdiff([1:length(obj.LAT.AclassinfoP)],inp.loc);
             LATnew=obj.LAT;
             try
             LATnew.AclabelP(locTrim)=[];
             end
             LATnew.AclassinfoP(locTrim)=[];
             LATnew.Apred(locTrim,:)=[];
             try
              LATnew.RawSpectra.Pset(locTrim,:)=[];   
             end
             out_obj=ssds(LATnew);
             
%              inp.pathfname_AT=textual_replaceBetween(obj.pathfname_AT,'_nsampP','_',num2str(length(inp.loc)));
             inp.pathfname_AT=strrep_keyword_between_markers_wlistRHS(obj.pathfname_AT,'_nsampP',{'_','.mat'},num2str(length(inp.loc)),'keepBoth');
             out_obj.saveAT(inp);
             
             
             
         end % end of extract_Pset method
        %%%%%%%%%%%%%%%%%%
        % generate new subtype of LAT by providing AclabelT_LC
        % create clistclslabel_LC and calculate AclassinfoMap2LC
        function out_obj=Label_Combination(obj,inp)
        % see also Label_Combination_Mapping and LC    
          disp('generate new subtype of LAT with "AclassinfoMap2LC"')    
            
        out_obj='';    
        end
        %%%%%%%%%%%%%%%%%%
        % generate new subtype of LAT by providing AclabelT_LC
        % create clistclslabel_LC and calculate AclassinfoMap2LC
        function out_obj=LC(obj,inp)
        % alias for Label_Combination
        % see also Label_Combination_Mapping and Label_Combination    
         error('still under construction')   
           out_obj=Label_Combination(obj,inp); 
            
        end
        %%%%%%%%%%%%%%%%%%
        function out_obj=asmc_Global_extract_Local(obj,inp)
            % apply auto based on Global (orig) classes
            % then extract local classes
            % see ppt: "SCiO OTC API ncls2 vs ILM vs ncls10_wPCA_SVM_autoGlobal_vs_autoLocal.pptx"
            % 
        %  need inp: e.g   inp.cls_pick={'C01','C03'};
        % see also: asmc_Orig_Tset_extract_subset() merge_rm_extract_class()
        % 
        inp.action='extract';
          out_obj=  asmc_Orig_Tset_extract_subset(obj.pathfname_AT,inp);
        
        end % end of asmc_Global_extract_Local()
        %%%%%%%%%%%%%%%%%%
        % split Atrainpk into odd vs even for cross validation etc operations
        function [out_obj_Todd_Peven out_obj_Teven_Podd]=Split_Odd_Even(obj)
        % after split, all Atrainpk samples used exactly once in all folds
        % modified from Atrainpk_Split_Odd_Even(pathfname_AT)
        % see also: test_ssds_Split_Odd_Even_method()
        pathfname_AT=obj.pathfname_AT;
        out=Atrainpk_Split_Odd_Even(pathfname_AT);
        
        out_obj_Todd_Peven=ssds(out.pathfname_Todd_Peven);
        out_obj_Teven_Podd=ssds(out.pathfname_Teven_Podd);

        
        end % end of Split_Odd_Even() method
        %%%%%%%%%%%%%%%%%%
        % split Atrainpk into Nfolds for cross validation etc operations
        function [array_obj_TP_pairs]=split_Nfolds(obj,inp)
        % after split, all Atrainpk samples used exactly once in all folds
        % modified from Atrainpk_Split_Odd_Even() and PLS_Tcv()
        % see also: test_ssds_split_Nfolds_method()
        error('under construction')
        end % end of split_Nfolds() method
        %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end  % end for methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end  % end for classdef