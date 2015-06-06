function varargout = syssli(varargin) %#ok<STOUT>
%%
z = {}; % NOTE: do not delete
%% ------ define local defined variables
gui_fig_color = [0.85 0.915 0.886];
gui_plotconfig_color = [0.792 0.906 0.976];
gui_importexport_color = [0.792 0.806 0.676];
gui_void_color = [0.9 0.9 0.9];

%--
legend_visibility = 'hide';             % show or hide
max_figplots = 30;                      % maximum seperate plot axis (default - 20)
num_param_sliders = 9;                  % number of slider paramters (currently can only choose 9)
invals = num_param_sliders;

h_slider = zeros(1,num_param_sliders);  %
slider_static = h_slider;  %
h_slider_edit_max = zeros(1,num_param_sliders);
h_slider_edit_min = zeros(1,num_param_sliders);
h_slider_edit_value = zeros(1,num_param_sliders);
h_slider_edit_name = zeros(1,num_param_sliders);
h_slider_button_reset = zeros(1,num_param_sliders);
h_slider_button_remove = zeros(1,num_param_sliders);
h_slider_button_removeandreset = zeros(1,num_param_sliders);


%-- axis events
IsMoveAxis = 0;
IsMoveAxis_orig_pos = [0 0];
IsMoveAxis_Xlim = [0 1];
IsMoveAxis_Ylim = [0 1];
XAxisFixed = 0; % default: false
YAxisFixed = 0; % default: false
zoomSpeed = 0.2;

%-- out save load
ExportImportParamDir = pwd;
ExportImportFigDir = pwd;
ExportImportGUIDir = pwd; %#ok<NASGU>


params = struct;                        % structure of parameter names with values
full_names = fieldnames(params);         % names of all changing parameters
% original_params = params;
params_static = params;                 % statically define a set of params incase of reset
ext = struct;

%--
emp_c = repmat({''},1,max_figplots);
func = struct('plot_set',{emp_c},'plot_command',{emp_c},'plot_title',{emp_c},'in','');
slider_names = repmat({''},1,num_param_sliders);    % names of slider  parameters
is_slider_active = zeros(1,num_param_sliders);         % values of corresponding parameter of slider


%slidermove = 0;


%% --- initialise the gui settings etc
    function initialise (handles) %#ok<INUSD>
        
        %slidermove = @(choice)slidermove_plot_textform(choice);

        if nargin < 1
            handles = struct; %#ok<NASGU>
        end
        if isempty (varargin)
            % if input is empty, this example will be produced
            % note: this is not an example to work from.
            % usser-defined variables to GUI variable are directly set here
            func.in = @example1;
            func.plot_set{1} = 'plot(z.range,z.sol,''DisplayName'',''test'');';
            func.plot_title{1} = 'sine wave';
            % set parameters to change
            params =  struct('x',10,'y',2,'jojo',1.0);
            % set names for slider initially and ammount of them (invals)
            slider_names{1}='x'; slider_names{2}='y'; invals = 2;
        else
            
            % the main file function
            func.in = varargin{1};
            
            % if arguments are input as structure [default]
            if isstruct(varargin{2})
                
                params = {};
                if isfield(varargin{2},'initial')
                    params = varargin{2}.initial;
                end
                
                full_names = fieldnames(params);
                
                if isfield(varargin{2},'init_slider') % if init slider exists
                    invals = numel(varargin{2}.init_slider); %nargin(fin)
                    if invals > num_param_sliders
                        warning(['init_slider argument can only have upto 9 parameter names, you have chosen ', numel(varargin{2}.init_slider), '. init_slider truncated to ', num_param_sliders, '.']);
                        invals = num_param_sliders;
                    end
                    slider_names(1:invals) = varargin{2}.init_slider(1:invals);
                else % ... if not: auto set min(params,9)
                    invals = min(num_param_sliders,length(full_names));
                    slider_names(1:invals) = full_names(1:invals);
                end
                
                % store external variables which you wish to store
                % [stored in .ext variable struct
                if isfield(varargin{2},'ext')
                    ext = varargin{2}.ext; %#ok<SETNU>
                end
                
                func.in = varargin{1};
                mStr = {'plot_title','plot_set'};
                if isfield(varargin{2},'plot_set')
                    % file version of plot list
                    if ischar(varargin{2}.plot_set)
                        parseFileToPlots (varargin{2}.plot_set)
                    else
                        for ch = 1:length(mStr)
                            if isfield(varargin{2},mStr{ch})
                                func.(mStr{ch})(1:numel(varargin{2}.(mStr{ch}))) =  varargin{2}.(mStr{ch});
                            else

                            end
                        end
                    end
                else
                    %error('no plot_set variable')
                end
                
            end
        end % end of isempty (varargin)
        %-
        for ch=1:num_param_sliders
            func.plot_command{ch} = displaystr_to_plotstr(strrep(func.plot_set{ch},';','; '));
            func.plot_set{ch} = plotstr_to_displaystr(strrep(func.plot_set{ch},';','; '));
        end
        
        %- organise slider-gui related properties
        for snum=1:invals
            set(h_slider_button_reset(snum), 'enable', 'on')
            set(h_slider_button_remove(snum), 'enable', 'on')
            set(h_slider_button_removeandreset(snum), 'enable', 'on')
            
            mval = params.(slider_names{snum});
            set(h_slider_edit_value(snum),'string',num2str(mval));
            set(h_slider_edit_max(snum),'string',num2str(max(mval*2,0) + (mval==0)*0.1));
            set(h_slider_edit_min(snum),'string',num2str(min(mval*2,0) - (mval==0)*0.1));
            
            adjust_sliderlevels(snum)
            set(h_slider_edit_name(snum),'string',slider_names{snum});
        
            % binary is slider active
            is_slider_active(snum) = 1; 
            % set static params
            %slider_static(snum) =  params.(slider_names{snum});
        end

        
        %try set(handles.edit_out,'string',func.out{1});
            %try set(handles.edit_out,'string',a,'position',b);
        %catch
        %    disp('the value for edit_out (i.e. plots) didnt work. you sure it was a string?')
        %end
        set(h_edit_out,'string',func.plot_set{1});
        set(h_edit_name,'string',func.plot_title{1});

        % post-initialisation - dependent on user input (above was standard
        % protocol if use did not specify certain arguments
        if invals < num_param_sliders - 1
            set(h_slider(invals + 1), 'enable','off'); %
        end
        if invals < num_param_sliders
            for k = invals+1:num_param_sliders
                set(h_slider_edit_value(k), 'enable','off')
                set(h_slider_edit_min(k), 'enable','off')
                set(h_slider_edit_max(k), 'enable','off')
                set(h_slider(k), 'enable','off')
                set(h_slider(k), 'BackgroundColor',[0.9 0.9 0.9])
            end
        end
        
        %- h_pop_plots - set names as stored titles
        tempA=get(h_pop_plots,'String');
        for iter=1:numel(func.plot_title)
            tempA{iter} = func.plot_title{iter};
        end
        set(h_pop_plots,'String',tempA);
        
        %- JAVA options for active horizontal scroll
        % source: http://undocumentedmatlab.com/blog/customizing-listbox-editbox-scrollbars
        jEdit = findjobj(h_edit_out);
        jEditbox = jEdit.getViewport().getComponent(0);
        jEditbox.setWrapping(false);                %# turn off word-wrapping
        set(jEdit,'HorizontalScrollBarPolicy',30);  %# HORIZONTAL_SCROLLBAR_AS_NEEDED
        % maintain horizontal scrollbar policy which reverts back on component resize
        hjEdit = handle(jEdit,'CallbackProperties');
        set(hjEdit, 'ComponentResizedCallback',...
            'set(gcbo,''HorizontalScrollBarPolicy'',30)')

        
        adjust_sliderlevels(1)
        set(hMain,'Units','Normalized')
        
        full_names = fieldnames(params);
        params_static = params;
                
%         k=whos
%         for niter=1:length(k)   
%             %k(niter).nesting
%             fprintf('%s : ',k(niter).name)
%             k(niter).nesting
%         end
    end


%% ------ parse text file into model for plotting
    function parseFileToPlots
        fid = fopen('dodo.txt');
        plotTxt = textscan(fid,'%s','Delimiter','\n');
        plotTxt=plotTxt{1};
        setIter=0;
        strFlush='';
        for lIter = 1:length(plotTxt)
            line = plotTxt(lIter);
            if strcmp(line(1:3),':::')
                setIter = setIter + 1;
                func.plot_title{setIter} = strtrim(line(4:end));
                if setIter > 1
                    func.plot_set{setIter-1} = strFlush;
                    strFlush = '';
                end
            else
                line = trim(line); % trim trailing and leading spaces
                if ~strcmp(line(end),';')
                    line(end+1) = ';'; %#ok<AGROW>
                end
                strFlush = strFlush + line;
            end
        end
    end
%% ------ defined changes to slider function during use
%% --
    function adjust_sliderlevels (choice)
        
        n_min = str2double(get(h_slider_edit_min(choice), 'string'));
        if ~isnumeric(n_min) || ~isfinite(n_min),
            n_min = 0;
        end
        n_max = str2double(get(h_slider_edit_max(choice), 'string'));
        if ~isnumeric(n_max) || ~isfinite(n_max),
            n_max = 0;
        end    
        n_val = str2double(get(h_slider_edit_value(choice), 'string'));
        if ~isnumeric(n_val) || ~isfinite(n_val),
            n_val = 0;
        end
        if n_val > n_max, n_max = n_val;
            set(h_slider_edit_max(choice),'string',num2str(n_max));
        elseif n_val < n_min, n_min = n_val;
            set(h_slider_edit_min(choice),'string',num2str(n_min));
        end
        set(h_slider(choice), 'min', n_min, 'max', n_max,'value',n_val)
        
    end


%-- FUNCTION: to add and remove hold on; hold off if not there
    function outstr = plotstr_to_displaystr (instr)
        % remove hold on                    - change first ';' to ; \nhold on;\n'
        instr = regexprep(instr,'hold on;','',1);
        % remove hold off                   - change last ';' to ; \nhold off;\n'
        instr = regexprep(instr,'hold off;','',1);
        % remove spaces
        % add new line                      - change ';' to ';\n'
        outstr = strrep( instr, ';', sprintf(';\n'));
    end
%-- FUNCTION: to add and remove hold on; hold off if not there
    function outstr = displaystr_to_plotstr (instr)
        % remove new lines              - remove '\n' to ''
        instr = strrep(sprintf(instr),sprintf('\n'),';');
        % add hold on                   - change first ';' to ; hold on;'
        instr = regexprep(instr,';',';hold on;',1);
        % add hold off                  - change last ';' to ; hold off;'
        outstr = regexprep(instr,';',';hold off;',numel(strfind(instr, ';')));
    end

%% ------ define gui window related functions : UPON USER COTROL
% --
%     function create_plot_edit(hObject, eventdata, handles)
%         if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%             set(hObject,'BackgroundColor','white');
%         end
%     end
%     function creation_edit_value_slider (hObject, eventdata, choice)
%         adj_sliderlevels(choice);
%         slidermove;    
%     end

%% ------ define gui window related functions : UPON CREATION


%-- slider move
    function slidermove (choice) %_plot_textform (choice)
        %persistent z
        p=get(h_slider,'value');
        
        if nargin<1; choice = 0; end
        %slider_names
        if choice > 0;
            params.(slider_names{choice}) = p{choice}(1);
        else
            for iter = find(is_slider_active==1)
                params.(slider_names{iter}) = p{iter}(1);
            end
        end
%         set(h_axes,'LegendColorbarListeners',[]); 
%         setappdata(h_axes,'LegendColorbarManualSpace',1);
%         setappdata(h_axes,'LegendColorbarReclaimSpace',1);
        %set(h_display_error,'String','');
        if ~isequal(slider_static,[p{:}]) || choice==-1 % || handles.profile == 1
        set(h_display_error,'String','');   
            if choice >= 0
                set(h_running_text,'visible','on')
                set(h_display_error,'String','');
                try  z = func.in(params,varargin{3:end});
                catch
                    e_str = sprintf('Error in @%1$s:\n%2$s',func2str(func.in),lasterr);
                    set(h_display_error,'String',e_str)
                end
                set(h_running_text,'visible','off');
            end
            xlim_statify = get(h_axes,'Xlim');
            ylim_statify = get(h_axes,'Ylim');
            h_out = func.plot_command{get(h_pop_plots,'value')};
            
            try 
                eval(h_out);
            catch
                %disp(lasterr);
                cla(h_axes);
                e_str = sprintf('Error in current plot configuration:\n%1$s',lasterr);
                set(h_display_error,'String',e_str)
            end
            if XAxisFixed == 1
                set(h_axes,'Xlim',xlim_statify);
            end
            if YAxisFixed == 1
                set(h_axes,'Ylim',ylim_statify);
            end
            %drawnow;
            
            for k = find(is_slider_active==1)
                set(h_slider_edit_value(k),'string',num2str(get(h_slider(k),'value')));
            end
            
        end
        %bb=tic
        if ~isempty(get(h_axes, 'children'))
            legend(h_axes,legend_visibility,'Location','best');
        end
        
        % update new static parameters upon success
        slider_static = [p{:}];
    end


%-- change plots from drop down menu
    function change_plots (varargin)
        % change title/desc of plot
        set(h_edit_name,'String',func.plot_title{get(h_pop_plots,'Value')})
        % change plot commands set
        set(h_edit_out,'String',func.plot_set{get(h_pop_plots,'Value')})
        if isempty(func.plot_set{get(h_pop_plots,'Value')})
            cla(h_axes,'reset') % clear
            % if this plot empty, send out warning
            set(h_display_error,'String','WARNING: the plot is empty, write and save or choose another plot set')
        else
            save_and_plot % reset clf
        end
    end

%-- update plot window for current
    function save_and_plot (varargin)
        % store new plotter title, if available
        func.plot_title{get(h_pop_plots,'Value')} = get(h_edit_name,'String');        
        func.plot_set{get(h_pop_plots,'Value')} = get(h_edit_out,'String');
        disp_str = get(h_edit_out,'String');
        plot_str = '';
        for iter=1:size(disp_str,1)
            plot_str = strcat(plot_str, disp_str(iter,:));
        end
        % change display version ready for plotting
        func.plot_command{get(h_pop_plots,'Value')} = displaystr_to_plotstr(plot_str);        
        if isempty(func.plot_command{get(h_pop_plots,'Value')})
            disp('boom')
            cla
        else
            slidermove( -1)
        end 
    end


%-- change the visibility of the legend
    function switch_legend_visibility (varargin)
        if get(h_check_legend,'value') == 0 % if not ticked
            legend_visibility='hide';
        else % if ticked (i.e. legend is active)
            legend_visibility='show';
        end
    end

%--
    function switch_checkbox_cursormode (varargin)
        if get(h_cursormode,'value')==0
            %handles.cursorMode.removeAllDataCursors()
            datacursormode(hMain,'off')
        else
            datacursormode(hMain,'on')
        end
    end

%-- movement of slider return callback to modify ...
    function active_slider_change (~, ~, choice)
        slidermove (choice);
    end

%-- edit minimum value of the respective parameter value
    function adjust_slider_limit(~, ~, choice)
        adjust_sliderlevels (choice)
    end

%-- edit maximum value of the respective parameter value
    function edit_name_slider(~, ~, choice, is_reset)
        if nargin < 4; is_reset = 0; end
        
        new_sname = get(h_slider_edit_name(choice),'string');
        orig_sname = slider_names{choice};
        
        %-
        if any(ismember(full_names,new_sname))~=0
            if strcmp(new_sname,orig_sname) && is_reset == 0
                set(h_display_error,'String',strcat('WARNING: The name ',new_sname,' already set here'))
            elseif any(ismember(slider_names,new_sname)) && is_reset == 0
                set(h_display_error,'String',strcat('WARNING: The name ',new_sname,''' already exists as dynamic parameter. Field changed back to ',orig_vname))
                set(h_slider_edit_name, 'string', orig_sname)
            else
                is_slider_active (choice) = 1;
                
                if is_reset == 1
                   params.(orig_sname) =  params_static.(orig_sname);
                end
                
                mval = params.(new_sname);
                set(h_slider_edit_value(choice), 'Enable','on','String', num2str(mval))
                set(h_slider_edit_max(choice), 'Enable','on','String', num2str(num2str(max(mval*2,0) + (mval==0)*0.1)))
                set(h_slider_edit_min(choice), 'Enable','on','String', num2str(num2str(min(mval*2,0) - (mval==0)*0.1)))
            
                slider_names{choice} = new_sname;
                set(h_slider_button_reset(choice), 'Enable', 'on')
                set(h_slider_button_remove(choice), 'Enable', 'on')
                set(h_slider_button_removeandreset(choice), 'Enable', 'on')

                adjust_sliderlevels(choice)
                set(h_slider(choice), 'BackgroundColor', [1 0.8 1], 'Enable', 'on')
            end
        elseif isempty(new_sname)
            set(h_slider_edit_name(choice),'String',slider_names{choice})
            set(h_display_error,'String',strcat('The string you supplied is empty. If you wish to remove this parameter, click the remove button'))            
        elseif any(ismember(full_names,new_sname))==0
            set(h_slider_edit_name(choice),'String',slider_names{choice})        
            set(h_display_error,'String',strcat('ERROR: The name ',new_sname,'''  does not exist in your parameter struct, Field changed back to ',orig_sname))
        end
        
    end

%-- empty slider cotainer
    function empty_slider (choice, is_reset)
        is_slider_active(choice) = 0;
        
        % set red slider for voidness
        set(h_slider(choice), 'BackgroundColor', gui_void_color, 'Enable', 'off')
        set(h_slider_edit_value(choice), 'String', '', 'Enable', 'off')
        set(h_slider_edit_min(choice), 'String', '0', 'Enable', 'off')
        set(h_slider_edit_max(choice), 'String', '0', 'Enable', 'off')      
    
        set(h_slider_button_reset(choice), 'Enable', 'off')
        set(h_slider_button_remove(choice), 'Enable', 'off')
        set(h_slider_button_removeandreset(choice), 'Enable', 'off')
    end


%-- statically edit the value of a particular slider variable
    function edit_value_slider (~, ~, choice)
        adjust_sliderlevels(choice);
        slidermove;    
    end

%-- remove and reset (to original value) the respective parameter slider
    function slider_param_removeandreset (hObject, eventdata, choice)
        edit_name_slider(hObject, eventdata, choice, 1)
        slider_param_remove (hObject, eventdata, choice)
    end

%-- remove the respective parameter slider
    function slider_param_remove (~, ~, choice)
        empty_slider(choice)
        set(h_slider_edit_name(choice), 'string', '')
        is_slider_active(choice) = 0;
        slider_names{choice} = '';
        
    end

%-- reset (to original value) the respective parameter slider
    function slider_param_reset (hObject, eventdata, choice)
    	edit_name_slider(hObject, eventdata, choice, 1)
        %adjust_sliderlevels (choice);
    end

%-- lock the x axis on graph from mouse move or plot update
    function lock_xaxis (varargin)
        XAxisFixed = get(h_fix_xaxis,'value');   
    end

%-- lock the y axis on graph from mouse move or plot update
    function lock_yaxis (varargin)
    	YAxisFixed = get(h_fix_yaxis,'Value');         
    end

%-- Change the slider zoom speed 
    function slider_zoom_change (varargin)
        zoomSpeed = 0.2 .* get(h_zoom_speed_slider,'Value');
    end

%-- save parameters
    function save_parameters (varargin)
        [file_islider, pathname_islider] = uiputfile(...
            {'*.m', 'program files (*.m)';...
            '*.mat','MAT-files (*.mat)'},...
            'Save current parameters as',...
            ExportImportParamDir);
        
        if file_islider ~= 0
            slidermove (0)
            pos_fs = strfind(file_islider,'.');
            file_islider_name = file_islider(1:pos_fs-1);
            file_islider_ext = file_islider(pos_fs+1:end);
            if strcmp(file_islider_ext,'m')
                fid_syssli = fopen(strcat(pathname_islider,file_islider),'w');
                fprintf(fid_syssli, 'function out = %s (out)\n', file_islider_name);
                fprintf(fid_syssli, 'if nargin == 0\n'); 
                fprintf(fid_syssli, '\tout = struct;\n'); 
                fprintf(fid_syssli, 'end\n\n'); 
                for c = 1:length(fieldnames(params))
                    fprintf(fid_syssli,'out.%s = %.16f;\n',full_names{c},params.(full_names{c}));
                end
                check=fclose(fid_syssli);
                clear fid_syssli
            elseif strcmp(file_islider_ext,'mat')
                save(strcat(pathname_islider,file_islider),'params')
            else
                set(handles.disperr,'String','SAVE PARAMETERS: File is either ''*.m'' or ''*.mat'' file. ')
            end
            disp('saved parameters')
        end
        
    end

%-- upload parameter
    function upload_parameters (varargin)
        [file_islider, pathname_islider] = uigetfile(...
            {'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Upload parameters',ExportImportParamDir);
        %ExportImportParamDir = pathname_islider;
        temp_params = params;
        if file_islider ~= 0
            file_islider_struct = strsplit(file_islider,'.');
            file_islider_name = file_islider_struct{1};
            file_islider_ext = file_islider_struct{2};

            if strcmp(file_islider_ext,'m')
                addpath(pathname_islider)
                eval(['params = ' file_islider_name ';']);
                rmpath(pathname_islider)
            elseif strcmp(file_islider_ext,'mat')
                loadmat_selection = matfile(strcat(pathname_islider,file_islider));
                params = loadmat_selection.params;
            end
            % check if 'params' variable is a struct params and compares to original
            % params
            %     A=fieldnames(orderfields(params));
            %     B=fieldnames(orderfields(orig_params));
            %     isequal(A,B)
            if isequal(fieldnames(orderfields(temp_params)),...
                    fieldnames(orderfields(params)))
                disp('uploaded file')
                for i=1:num_param_sliders % was invals
                    if get(h_slider(i),'enable','on')
                        mval = params.(slider_names{i});
                        set(h_slider_edit_value(i),'string',num2str(num2str(mval)));
                        set(h_slider_edit_max(i),'string',num2str(max(mval*2,0) + (mval==0)*0.1));
                        set(h_slider_edit_min(i),'string',num2str(min(mval*2,0) - (mval==0)*0.1));
                        
                        adjust_sliderlevels(i);
                    end
                end
                slidermove;
                params_static = params;
            else
                % Since the two param structures dont have the same fieldnames, an
                % error is displayed. In this version it is not possible to overlap
                % 'params' file. In future releases, this may be possible
                set(h_display_error,'String','ERROR in UPLOAD PARAMS: param file does not match with current open file')
                % is_paramfile_struct = isstruct(params);
            end
        else
            % not supposed to reach this
        end        
    end

%-- upload parameter
    function save_figure (varargin)
        [file_islider, pathname_islider] = uiputfile(...
            {'*.mat','MAT-files (*.mat)';
            '*.fig','FIG-files (*.fig)'},...
            'Save figure as',...
            ExportImportFigDir);
        %ExportImportFigDir = pathname_islider;

        if file_islider ~= 0
            file_ext_is_last_word = strsplit(file_islider,'.');
            if strcmp(file_ext_is_last_word(end),'fig')
                Fig_temp = figure('Visible','on');
                Axes_temp = h_axes;
                copyobj(Axes_temp, Fig_temp);
                findall(Fig_temp,'type','axes')
                set(findall(Fig_temp,'type','axes'),'Units','Normalized','position',get(0,'factoryAxesPosition'))
                savefig(Fig_temp,strcat(pathname_islider,file_islider));
                close(Fig_temp);
                disp('Figure axis data saved [.fig]...')
            elseif strcmp(file_ext_is_last_word(end),'mat')
                Fig_temp = figure('Visible','on');
                Axes_temp = h_axes;
                copyobj(Axes_temp, Fig_temp);
                findall(Fig_temp,'type','axes')
                set(findall(Fig_temp,'Units','Normalized','type','axes'),'position',get(0,'factoryAxesPosition'))
                savefig(Fig_temp,'_temporary_figurex.fig');
                close(Fig_temp);
                
                varstore = load('_temporary_figurex.fig', '-mat'); %#ok<NASGU>
                delete('_temporary_figurex.fig')
                save(strcat(pathname_islider,file_islider),'varstore')
                disp('Figure axis data saved [.mat]...')
            end
        end
        
    end

%-- save gui state
    function save_GUI_state (varargin)
        %h_save_gui_state
        [file_islider, pathname_islider] = uiputfile(...
            {'*.sli','SYSLIDER state (*.sli)'},...
            'Save current state of SYSSLIDER as...');
        if file_islider ~= 0
%             jo=whos;
%             sliderprops(1).nesting.function
%             jo(1).nesting.level
%             jo
%             %type(jo.nesting)
%             mfilename
%             
%             for c = 1:length(jo)
%                 c
%             end
            
%            pos_fs = strfind(file_islider,'.');
%            f_ext = file_islider(pos_fs+1:end);
%             if ~strcmp(ExportImportGUIDir,pathname_islider(1:end-1))
%                 set(h_display_error,'string','ERROR in GUI state save, can only save in current directory from which sysslider is run')
%             elseif strcmp(f_ext,'fig')
%                 disp('GUI state saved...')
%                 savefig(hMain,strcat(pathname_islider,file_islider));
%             else
%             end
        end
    end




%% ------ define GUI window objects
%% --- define figure window properties
hMain = figure(...
    'Units','Pixels',...
    'PaperUnits',get(0,'defaultfigurePaperUnits'),...
    'IntegerHandle','off',...
    'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
    'Name','Computational & Systems Biophysics : SYSLI',...
    'NumberTitle','off',...
    'PaperPosition',get(0,'defaultfigurePaperPosition'),...
    'PaperSize',get(0,'defaultfigurePaperSize'),...
    'PaperType',get(0,'defaultfigurePaperType'),...
    'Position',[38 66 950 730],...
    'HandleVisibility','callback',...
    'Tag','sysslider_fig',...
    'UserData',[],...
    'Visible','on',...
    'resize','off',...
    'WindowButtonMotionFcn',@mouse_motion_update,...
    'WindowButtonDownFcn',@mouse_button_activate_axis_move,...
    'WindowButtonUpFcn',@mouse_button_deactivate_axis_move,...
    'WindowScrollWheelFcn',@mouse_scroll_zoom_axis,...
    'CloseRequestFcn',@window_closeFcn,...
    'color',gui_fig_color);
%     'ToolBar','none',...
%     'MenuBar','none',...
    
%% -- define slider configeration settings
properties_hMain_text = {...
    'Parent',hMain,...
    'Units','normalized',...
    'FontName',get(0,'defaultuicontrolFontName'),...
    'FontSize',get(0,'defaultuicontrolFontSize'),...
    'Style','text',...
    'backgroundcolor',get(hMain,'color')};

h_running_text = annotation(...
    hMain,...
    'textbox',[0.17 0.684 0.263 0.1],...
    'Color',[0.7 0.3 0.3],...
    'String','Loading...',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'visible','off',...
    'FontSize',45.0);
% h_running_text = uicontrol(...
%     properties_hMain_text{:},...
%     'Position',[0.19 0.685 0.255 0.1],...
%     'String','Loading...',...
%     'FontSize',40.0,...
%     'ForegroundColor','red',...
%     'BackgroundColor','white',...
%     'Tag','text_loading');

uicontrol(...
    properties_hMain_text{:},...
    'Position',[0.775 0.425 0.063 0.02],...
    'String','Min',...
    'Tag','text5');

uicontrol(...
    properties_hMain_text{:},...
    'Position',[0.85 0.425 0.063 0.02],...
    'String','Max',...
    'Tag','text6');

uicontrol(...
    properties_hMain_text{:},...
    'Position',[0.925 0.425 0.063 0.02],...
    'String','Value',...
    'Tag','text7');

properties_slider_universal = {...
    'Parent',hMain,...
    'Units','normalized',...
	'FontName',get(0,'defaultuicontrolFontName'),...
	'FontSize',get(0,'defaultuicontrolFontSize')};

for schoice = 1:num_param_sliders
    h_slider(schoice) =  uicontrol(...
        properties_slider_universal{:},...
        'BackgroundColor',[1 0.8 1],...
        'Callback',@(hObject,eventdata)active_slider_change(hObject,eventdata,schoice),...% slider_Callback(schoice),...
        'CData',[],...
        'Position',[0.09 0.39-0.034*schoice+0.034 0.53 0.032],...
        'String',{  'ƒXƒ‰ƒCƒ_' },...
        'Style','slider',...
        'Tag',strcat('slider',num2str(schoice)),...
        'UserData',[]);
    %        'KeyPressFcn',{@local_CreateFcn, @(hObject,eventdata)sysslider(strcat('slider_KeyPressFcn'),hObject,eventdata,guidata(hObject)), appdata},...%slider_KeyPressFcn
    
    h_slider_edit_min(schoice) = uicontrol(...
        properties_slider_universal{:},...
        'BackgroundColor',[0.7 1 0.7],...
        'Callback',@(hObject,eventdata)adjust_slider_limit(hObject,eventdata,schoice),...%edit_min_Callback
        'Position',[0.775 0.39-0.034*schoice+0.034 0.065 0.034],...
        'String','0',...
        'Style','edit',...
        'Tag',strcat('edit_min',int2str(schoice)));
    
    h_slider_edit_max(schoice) = uicontrol(...
        properties_slider_universal{:},...
        'BackgroundColor',[0.7 1 0.7],...
        'Callback',@(hObject,eventdata)adjust_slider_limit(hObject,eventdata,schoice),...%edit_min_Callback
        'Position',[0.85 0.39-0.034*schoice+0.034 0.065 0.034],...
        'String','1',...
        'Style','edit',...
        'Tag',strcat('edit_max',int2str(schoice)));
    
    h_slider_edit_name(schoice) = uicontrol(...
        properties_slider_universal{:},...
        'BackgroundColor',[0.7 0.7 1],...
        'FontWeight','bold',...
        'Callback',@(hObject,eventdata)edit_name_slider(hObject,eventdata,schoice),...%edit_vname_Callback'
        'Position',[0.015 0.39-0.034*schoice+0.034 0.065 0.034],...
        'String','',...%strcat('x',int2str(schoice)),...
        'Style','edit',...
        'Tag',strcat('text_x',int2str(schoice)));
    
    h_slider_edit_value(schoice) = uicontrol(...
        properties_slider_universal{:},...
        'BackgroundColor',[1 0.7 0.7],...
        'Callback',@(hObject,eventdata)edit_value_slider(hObject,eventdata,schoice),... % edit_val_Callback
        'Position',[0.925 0.39-0.034*schoice+0.034 0.065 0.034],...
        'String',blanks(0),...
        'Style','edit',...
        'Tag',strcat('edit_val',int2str(schoice)));
    
    h_slider_button_reset(schoice) = uicontrol(...
        properties_slider_universal{:},...
        'Callback',@(hObject,eventdata)slider_param_reset(hObject,eventdata,schoice),... %edit_pushbutton_reset_Callback
        'Position',[0.63 0.39-0.034*schoice+0.034 0.04 0.032],...
        'String','Reset',...
        'enable','off',...
        'Tag',strcat('pushbutton_reset',int2str(schoice)));
    
    h_slider_button_removeandreset(schoice) = uicontrol(...
        properties_slider_universal{:},...
        'Callback',@(hObject,eventdata)slider_param_removeandreset(hObject,eventdata,schoice),... %edit_pushbutton_removeandreset_Callback
        'Position',[0.73 0.39-0.034*schoice+0.034 0.04 0.032],...
        'String','Both',...
        'enable','off',...
        'Tag',strcat('pushbutton_removeandreset',int2str(schoice)));
    
    h_slider_button_remove(schoice) = uicontrol(...
        properties_slider_universal{:},...
        'Callback',@(hObject,eventdata)slider_param_remove(hObject,eventdata,schoice),... %edit_pushbutton_remove_Callback
        'Position',[0.675 0.39-0.034*schoice+0.034 0.05 0.032],...
        'String','discard',...
        'enable','off',...
        'Tag',strcat('pushbutton_remove',int2str(schoice)));
end

%% --- define plot configeration section
hPan = uipanel(hMain, 'Title','Plot configurations', ...
    'Units','normalized', 'Position',[0.64 0.64 0.35 0.33],...
    'backgroundcolor',gui_plotconfig_color,'Tag','hPan');

properties_hPan_config = {...
    'Parent',hPan,...
    'Units','normalized',...
    'FontName',get(0,'defaultuicontrolFontName'),...
    'FontSize',get(0,'defaultuicontrolFontSize'),...
    };

h_edit_out = uicontrol(...
    properties_hPan_config{:},...
    'BackgroundColor',[1 1 1],...
    'FontSize',12,...
    'HorizontalAlignment','left', ...
    'Position',[0.02 0.12 0.94 0.55],...
    'String',blanks(0),...
    'Style','edit',...
    'max',10.0,...
    'Tag','edit_out');

h_edit_name = uicontrol(...
    properties_hPan_config{:},...
    'BackgroundColor',[1 1 1],...
    'FontSize',12,...
    'FontWeight','bold',...
    'ForeGroundColor',[0 0 1],...
    'HorizontalAlignment','left', ...
    'Position',[0.02 0.7 0.94 0.16],...
    'String',blanks(0),...
    'Style','edit',...
    'max',2.0,...
    'Tag','edit_name');

temp_plot_num_str = cell(1,max_figplots); 
temp_plot_num_str(:) = {'plots'};
for ipl=1:max_figplots
    temp_plot_num_str{ipl} = [temp_plot_num_str{ipl},' ',num2str(ipl)];
end
h_pop_plots = uicontrol(...
    properties_hPan_config{:},...
    'Position',[0.45 0.88 0.5 0.1],...
    'String',temp_plot_num_str,...
    'Style','popup',...
    'Tag','pop_plots',...
    'Callback',@change_plots,...
    'CreateFcn', blanks(0));

h_save_and_plot = uicontrol(...
    properties_hPan_config{:},...
    'Callback',@save_and_plot,...
    'Position',[0.02 0.89 0.4 0.1],...
    'String','save & plot',...
    'Tag','save_and_plot'); %#ok<NASGU>

h_check_legend = uicontrol(...
    properties_hPan_config{:},...
    'Callback',@switch_legend_visibility,...
    'Position',[0.02 0.02 0.2 0.07],...
    'String','Legend',...
    'Style','checkbox',...
    'Value',0.0,...
    'Tag','checkbox_legend',...
    'backgroundcolor',get(hPan,'backgroundcolor'));

h_cursormode = uicontrol(...
    properties_hPan_config{:},...
    'Callback',@switch_checkbox_cursormode,...
    'Position',[0.2 0.02 0.3 0.07],...
    'String','cursor',...
    'Style','checkbox',...
    'Value',0.0,...
    'Tag','checkbox_cursormode',...
    'backgroundcolor',get(hPan,'backgroundcolor'));

h_fix_xaxis = uicontrol(...
    properties_hPan_config{:},...
    'Callback',@lock_xaxis,... %edit_checkbox_XAxisFix
    'Position',[0.37 0.02 0.3 0.07],...
    'String','x-axis lock',...
    'Style','checkbox',...
    'Value',0.0,...
    'Tag','checkbox_XAxisFix',...
    'backgroundcolor',get(hPan,'backgroundcolor'));

h_fix_yaxis = uicontrol(...
    properties_hPan_config{:},...
    'Callback',@lock_yaxis,...%edit_checkbox_YAxisFix
    'Position',[0.60 0.02 0.3 0.07],...
    'String','y-axis lock',...
    'Style','checkbox',...
    'Value',0.0,...
    'Tag','checkbox_YAxisFix',...
    'backgroundcolor',get(hPan,'backgroundcolor'));

%% ----
h_axes = axes(...
    'Parent',hMain,...
    'Position',[0.04 0.50 0.54 0.46],...
    'CameraPosition',[0.5 0.5 9.16025403784439],...
    'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
    'Color',get(0,'defaultaxesColor'),...
    'ColorOrder',get(0,'defaultaxesColorOrder'),...
    'FontName',get(0,'defaultaxesFontName'),...
    'FontSize',get(0,'defaultaxesFontSize'),...
    'LooseInset',[0.0567146509183546 0.0398395126424734 0.0414453218249514 0.0271633040744137],...
    'XColor',get(0,'defaultaxesXColor'),...
    'YColor',get(0,'defaultaxesYColor'),...
    'ZColor',get(0,'defaultaxesZColor'),...
    'Tag','axes1');


%% ---- export / import section
hPanUpDownLoad = uipanel(hMain, 'Title','Export/import parameters', ...
    'Units','normalized', 'Position',[0.64 0.005 0.35 0.1],...
    'backgroundcolor',gui_importexport_color);

properties_importexport_config = {...
    'Parent',hPanUpDownLoad,...
    'Units','normalized',...
    'FontName',get(0,'defaultuicontrolFontName'),...
    'FontSize',get(0,'defaultuicontrolFontSize'),...
    };

h_save_params = uicontrol(...
    properties_importexport_config{:},...
    'Callback',@save_parameters,...'edit_save_parameters',hObject,eventdata,guidata(hObject)),...
    'Position',[0.55 0.6 0.35 0.4],...
    'String','Save parameters'); %#ok<NASGU>

h_upload_params = uicontrol(...
    properties_importexport_config{:},...
    'Callback',@upload_parameters,...'edit_upload_parameters',hObject,eventdata,guidata(hObject)),...
    'Position',[0.55 0.1 0.35 0.4],...
    'String','Upload parameters',...
    'Tag','upload_params'); %#ok<NASGU>

h_save_fig = uicontrol(...
    properties_importexport_config{:},...
    'Callback',@save_figure,...'edit_save_figure',hObject,eventdata,guidata(hObject)),...
    'Position',[0.15 0.6 0.35 0.4],...
    'String','Save figure axis',...
    'Tag','save_figure'); %#ok<NASGU>

h_save_gui_state = uicontrol(...
    properties_importexport_config{:},...
    'Callback',@(hObject,eventdata)save_GUI_state(hObject,eventdata,guidata(hObject)),...%@save_GUI_state,...'edit_save_gui_state',hObject,eventdata,guidata(hObject)),...
    'Position',[0.15 0.1 0.35 0.4],...
    'String','Save GUI state',...
    'Tag','save_gui_state'); %#ok<NASGU>

%% -- warning.error dialog
appdata = [];
appdata.lastValidTag = 'hPanError';
hPanError = uipanel(hMain, 'Title','Warning/Error messages', ...
    'Units','normalized', 'Position',[0.64 0.50 0.35 0.12],...
    'backgroundcolor',[1 0.8 0.8]);

h_display_error = uicontrol(...
    'Parent',hPanError,...
    'Units','normalized',...
    'FontName',get(0,'defaultuicontrolFontName'),...
    'FontSize',get(0,'defaultuicontrolFontSize'),...
    'String','warning/error codes presented here',...
    'Style','text',...
    'Position',[0.02 0.02 0.96 0.96],...
    'Tag','disperr',...
    'horizontalalignment','left',...
    'max',2,...
    'fontsize',get(0,'defaultuicontrolFontSize')+2,...
    'ForegroundColor','red',...
    'BackgroundColor',get(hPanError,'backgroundcolor'));


%% --- additional unorganised features
h_zoom_speed_slider =  uicontrol(...
    'Parent',hMain,...
    'Units','normalized',...
	'FontName',get(0,'defaultuicontrolFontName'),...
	'FontSize',get(0,'defaultuicontrolFontSize'),...
    'BackgroundColor',[1 0.6 0.7],...
    'Callback',@slider_zoom_change,...
    'Position',[0.02 0.04 0.2 0.022],...
    'Style','slider',...
    'max',2.0,...
    'min',0.001,...
    'Value',1.0,...
    'Tag','sliderzoom');

%% --- exmaple function for testing purposes
    function out = example1(mParam)
        out.range = 1:0.1:100;
        out.sol = mParam.x*sin(mParam.y.*out.range);
    end


%--- windows figure mouse controls
    function mouse_motion_update (varargin)
        slidermove;
        if IsMoveAxis == 1
            mouse_current_point = get(varargin{1},'CurrentPoint');
            mouse_position = get(h_axes,'Position');
            % get ratio position moved on x axis relative to x axis position
            if XAxisFixed==0
                Xrange = IsMoveAxis_Xlim(2)-IsMoveAxis_Xlim(1);
                ratioX_posmove_axesrel = (IsMoveAxis_orig_pos(1) - mouse_current_point(1))/(mouse_position(3));
                set(h_axes,'XLim',IsMoveAxis_Xlim + ratioX_posmove_axesrel * Xrange);
            end
            % get ratio position moved on y axis
            if YAxisFixed==0
                Yrange = IsMoveAxis_Ylim(2)-IsMoveAxis_Ylim(1);
                ratioY_posmove_axesrel = (IsMoveAxis_orig_pos(2) - mouse_current_point(2))/(mouse_position(4));
                set(h_axes,'YLim',IsMoveAxis_Ylim + ratioY_posmove_axesrel * Yrange);
            end
        end
    end

%--
    function mouse_button_activate_axis_move (varargin)
        mouse_selectionType = get(varargin{1}, 'selectionType');
        mouse_current_point = get(varargin{1},'CurrentPoint');
        mouse_position = get(h_axes,'Position');
        if (mouse_current_point(1)>=mouse_position(1) && mouse_current_point(1)<= mouse_position(1)+mouse_position(3)&&...
        	mouse_current_point(2)>=mouse_position(2) && mouse_current_point(2)<= mouse_position(2)+mouse_position(4)&&...
            strcmp(mouse_selectionType, 'alt'))
                IsMoveAxis = 1; % is moveable
                IsMoveAxis_orig_pos = mouse_current_point;
                IsMoveAxis_Xlim = get(h_axes,'XLim');
                IsMoveAxis_Ylim = get(h_axes,'YLim');
%         elseif strcmp(mouse_selectionType, 'normal')
%             legend_visibility = 'hide'
        end
    end

%--
    function mouse_button_deactivate_axis_move (varargin)
        mouse_selectionType = get(varargin{1}, 'selectionType');
        if (strcmp(mouse_selectionType, 'alt'))
            IsMoveAxis = 0; % is moveable
%         elseif (strcmp(mouse_selectionType, 'normal'))
%             if get(h_check_legend,'value') == 1 % set not to hidden
%                 legend_visibility = 'show'
%             end
        end
    end

%--
    function mouse_scroll_zoom_axis (varargin)
        upOrDwn = varargin{2}.VerticalScrollCount; % from eventData
        mouse_current_point = get(varargin{1},'CurrentPoint');
        mouse_position = get(h_axes,'Position');
        X=mouse_current_point(1); Y=mouse_current_point(2);
        curXlim = get(h_axes,'XLim');
        curYlim = get(h_axes,'YLim');
        curXRange = curXlim(2)-curXlim(1);
        curYRange = curYlim(2)-curYlim(1);
        %if X>=0.04 && X<=0.62 && Y>= 0.47 && Y<=0.97
        if  X>=mouse_position(1) && X<=mouse_position(1)+mouse_position(3) && Y>=mouse_position(2) && Y<=mouse_position(2)+mouse_position(4)
            ratio_lowerX = (X-mouse_position(1))/mouse_position(3);
            ratio_upperX = 1-ratio_lowerX;
            ratio_lowerY = (Y-mouse_position(2))/mouse_position(4);
            ratio_upperY = 1-ratio_lowerY;
            if upOrDwn == -1
                %disp('up')
                if XAxisFixed==0
                    set(h_axes,'XLim',[curXlim(1)-zoomSpeed*curXRange*ratio_lowerX curXlim(2)+zoomSpeed*curXRange*ratio_upperX])
                end
                if YAxisFixed==0
                    set(h_axes,'YLim',[curYlim(1)-zoomSpeed*curYRange*ratio_lowerY curYlim(2)+zoomSpeed*curYRange*ratio_upperY])
                end
            else
                %disp('down')
                if XAxisFixed==0
                    set(h_axes,'XLim',[curXlim(1)+zoomSpeed*curXRange*ratio_lowerX curXlim(2)-zoomSpeed*curXRange*ratio_upperX])
                end
                if YAxisFixed==0
                    set(h_axes,'YLim',[curYlim(1)+zoomSpeed*curYRange*ratio_lowerY curYlim(2)-zoomSpeed*curYRange*ratio_upperY])
                end
            end
        end    
    end

%-- Window close callback clear variables
    function window_closeFcn (varargin)
        disp('syssli closed...')
        clear h_slider slider_static h_slider_edit_max = zeros(1,num_param_sliders);
        clear h_slider_edit_min h_slider_edit_value clear h_slider_edit_name
        clear params params_static full_names emp_c
        clear static_params func slider_names is_slider_active
        delete(get(0,'CurrentFigure'))
    end

    initialise;
end % end of sysslider function