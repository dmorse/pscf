function polymer_visual(filename)

close all
%% Reading the file
tic
C = textread(filename, '%s','delimiter', '\n');

ic = 0;
ndata = 0;
required_rep = 5;


while ndata <= required_rep
    ic = ic +1;
    
    if round(sum (sscanf(C{ic},'%f')),2)== 1.00 && round(sum (sscanf(C{ic+1},'%f')),2)== 1.00
        ndata = ndata+1;
    else
        ndata = 0;
    end
    
    if strcmp(strrep(char(C(ic)),' ', ''),'dim')==1
        dim = str2double(C{ic+1});          % Reads the grid dimensions
    elseif strcmp(strrep(char(C(ic)),' ', ''),'crystal_system')==1
        type = strrep(C(ic+1), '''', '');   % Reads the system type
    elseif strcmp(strrep(char(C(ic)),' ', ''),'cell_param')==1
        param = sscanf(C{ic+1},'%f')';            % Reads the cell parameters
    elseif strcmp(strrep(char(C(ic)),' ', ''),'N_monomer')==1
        n_mnr = str2double(C{ic+1});        % Reads the number of monomers
    elseif strcmp(strrep(char(C(ic)),' ', ''),'ngrid')==1
        grid = sscanf(C{ic+1},'%f')';            % Reads the grid size
    end
end

end_info = ic - required_rep - 1;                % Records the row in which the supplementary information ends
start_row = ic - required_rep;                  % Records the row in which the volume fractions start


%% Reading the grid points from the file

A = zeros(length(C) - end_info,n_mnr);

for i =start_row:length(C)
    A(i - end_info,:) = sscanf(C{i},'%f')';
end


%% Other Inputs


ncolour = 8; %Number of stored colourmaps (7 maps + 1 fill)
n_dp = 3; %Number of significant decimal places for colourmapping

inputvec = [1 1 1];
h_set = 0:3;
k_set = 0:1;
l_set = 0:1;
isovalue = ones(1,n_mnr)*0.25;
opacity = ones(n_mnr,2);
%opacity = [1,1;0,0.65;1,1];
%isovalue = [0.15,0.35,0.35];
isovalue = 'auto';
thick = 1; %Box Thickness Value
box_clr = [0.5 0.5 0.5]; %Box Colour
for in = 1:n_mnr
    mono_disp(in) = in;
    map_choice(in)= in;
    comp_disp(in) = in;
    weight(in) = 1;
end

%weight(1) = 1.2;
%weight(2) = 1.2;
%weight(3) = 1.2;
drawscatter = [1];
 %map_choice = [7 5 6]; %Map colours
%mono_disp = [];    %Desired monomers to display
%comp_disp = []; %Desired monomers to display in the composite



% C2 = textread('outputfileinput9.txt', '%s','delimiter', '\n');
%
% ic = 0;
%
% for ic = 1:length(C2)
%
%     if strcmp(strrep(char(C2(ic)),' ', ''),'thickness')==1
%         thick = str2double(C2{ic+1});          % Reads the box thickness
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'box_clr')==1
%         box_clr = sscanf(C2{ic+1},'%f');           % Reads the box colour
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'opacity')==1
%         opacity(:,1) = sscanf(C2{ic+1},'%f')'; % Reads the opacity
%         opacity(:,2) = sscanf(C2{ic+2},'%f')'; % Reads the opacity
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'isovalue')==1
%         if ischar(strrep(char(C2(ic+1)),' ', ''))
%         isovalue = strrep(char(C2(ic+1)),' ', '');
%         else
%         isovalue = sscanf(C2{ic+1},'%f')';           % Reads the isovalues
%         end
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'map_choice')==1
%         map_choice = sscanf(C2{ic+1},'%f')';    % Reads the chosen map colours
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'mono_disp')==1
%         mono_disp = sscanf(C2{ic+1},'%f')';           % Reads the desired monomers to display
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'comp_disp')==1
%         comp_disp = sscanf(C2{ic+1},'%f')';          % Reads the desired monomers to display in the composite
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'weight')==1
%         weight = sscanf(C2{ic+1},'%f')';          % Reads the desired weight of composite monomers
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'scatterdraw')==1
%         drawscatter = sscanf(C2{ic+1},'%f')';     % For scattering
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'h_set')==1
%         h_set = sscanf(C2{ic+1},'%f')';     % For scattering
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'k_set')==1
%         k_set = sscanf(C2{ic+1},'%f')';     % For scattering
%     elseif strcmp(strrep(char(C2(ic)),' ', ''),'l_set')==1
%         l_set = sscanf(C2{ic+1},'%f')';     % For scattering
%
%     end
% end



if strcmp(isovalue,'auto')==1
    linedraw = 1;
else
    linedraw=0;
end

%% Separating grid dimensions and angles

if strcmp(type,'hexagonal') == 1
    angle = [pi/2 pi/2 (2*pi)/3];
    cell_d = param;
elseif strcmp(type,'cubic') == 1
    angle = [pi/2 pi/2 pi/2];
    cell_d = param;
elseif strcmp(type,'tetragonal') == 1
    angle = [pi/2 pi/2 pi/2];
    cell_d = param;
elseif strcmp(type,'orthorhombic') == 1
    angle = [pi/2 pi/2 pi/2];
    cell_d = param;
elseif strcmp(type,'triclinic') == 1
    angle = [param(4) param(5) param(6)];
    cell_d = [param(1) param(2) param(3)];
elseif strcmp(type,'monoclinic') == 1
    angle = [pi/2 param(4) pi/2];
    cell_d = [param(1) param(2) param(3)];
elseif strcmp(type,'trigonal') == 1
    angle = [param(2) param(2) param(2)];
    cell_d = [param(1)];
elseif strcmp(type,'lamellar') == 1
    angle = [pi/2 pi/2 pi/2];
    cell_d = param;
else
    angle = [pi/2 pi/2 pi/2];
    cell_d = param;
end

%% Calculating the cell dimensions

if(length(cell_d)==1)
    new_cell = ones(1,3)*cell_d;          % Cubic crystals
elseif(length(param)==2)
    new_cell(1:2) = cell_d(1);            % Tetragonal crystals
    new_cell(3)   = cell_d(2);
else
    new_cell = cell_d;                    % Orthorhombic crystals
end

clear cell_d;cell_d = new_cell;

if(length(grid)==1)
    grid(2) = grid(1);                   % 3D grid for 1D crystals
    grid(3) = grid(1);
elseif(length(grid)==2)
    grid(3) = grid(1);                  % 3D grid for 2D crystals
end

for ig = 1:3
    plot_grid(ig) = grid(ig)+1;
end

%% Formulating the volume fraction in 4D  arrays for 3D visualization.

x = zeros(grid);
y = zeros(grid);
z = zeros(grid);
R = zeros([grid n_mnr]);
counter = 0;


for iz=1:grid(3)+1,
    for iy=1:grid(2)+1,
        for ix=1:grid(1)+1,
            counter = counter + 1;
            x(ix,iy,iz) = cell_d(1) * (ix-1)/grid(1) + (cos(angle(3)))*( cell_d(2) * (iy-1)/grid(2)) + ((iz-1)/grid(3))*(cos(angle(1))*cell_d(3));
            y(ix,iy,iz) = cell_d(2) * (iy-1)/grid(2) * sin(angle(3)) + ((iz-1)/grid(3))*cos(angle(2))*cell_d(3);
            z(ix,iy,iz) = cell_d(3) * (iz-1)/grid(3) * sin(angle(1)) * sin(angle(2));
            for in = 1:n_mnr
                if ix == grid(1)+1
                    R(grid(1)+1,:,:,in) = R(1,:,:,in);
                    counter = counter - (1/n_mnr);
                elseif iy == grid(2) + 1
                    R(:,grid(2)+1,:,in) = R(:,1,:,in);
                    counter = counter - (1/n_mnr);
                elseif iz == grid(3) + 1
                    R(:,:,grid(3)+1,in) = R(:,:,1,in);
                    counter = counter - (1/n_mnr);
                else
                    R(ix,iy,iz,in) = A(round(counter),in);
                end
            end
        end
        if(dim==1)
            counter = 0;
        end
    end
    if(dim==2)
        counter=0;
    end
end


%% Isovalue calculation

if strcmp(isovalue,'auto')==1
    
    isovalue = zeros(1,n_mnr);
    %Finding the range and limits
    n_lines = 1;
    
    for in = 1:n_mnr
        polmaxa(in) = max(max(max((R(:,:,:,in)))));
        polmin(in) = min(min(min((R(:,:,:,in)))));
        l_length(in) = polmaxa(in) -polmin(in);
    end
    
    
    %Creating the scaled matrix, S
    S = R;
    for in = 1:n_mnr
        for iz=1:grid(3)+1,
            for iy=1:grid(2)+1,
                for ix=1:grid(1)+1,
                    S (ix,iy,iz,in) = weight(in)*(R(ix,iy,iz,in)-polmin(in))/l_length(in); % Matrix with scaled values
                end
            end
        end
    end
    
    % Finding the location of the first (three) polmax points
    pmax_loc = zeros(n_lines,3,n_mnr);
    for in =1:n_mnr
        [si,sj,sk] = ind2sub([plot_grid n_mnr],find(S(:,:,:,in) == weight(in),n_lines));
        pmax_loc(:,:,in) = [si,sj,sk];
    end
    
    % Interleave
    point_series = zeros(n_lines*n_mnr,3);
    for ir = 1:n_lines
        for in = 1:n_mnr
            point_series(n_mnr*(ir-1)+in,:) = pmax_loc(ir,:,in); %Series of points to be plotted
        end
    end
    
    % Creating the lines
    clear line line_new x_plot intervalue ix iy iz
    l_size = 0;
    for ir = 1:(n_lines*n_mnr)-1 %Number of lines
        start_coord = point_series(ir,:);
        end_coord = point_series(ir+1,:);
        dir_vec(ir,:) = end_coord-start_coord;
        step_length(ir) = max(abs(dir_vec(ir,:)));
        for il = 1:step_length(ir)
            ix(ir,il) = start_coord(1)+ round((il-1)*(dir_vec(ir,1)/step_length(ir)));
            iy(ir,il) = start_coord(2)+ round((il-1)*(dir_vec(ir,2)/step_length(ir)));
            iz(ir,il) = start_coord(3)+ round((il-1)*(dir_vec(ir,3)/step_length(ir)));
            x_plot(il+l_size) = il+l_size;
            for in= 1:n_mnr
                line (in,il+l_size) = R (ix(ir,il),iy(ir,il),iz(ir,il),in);
                
            end
        end
        l_size = l_size + step_length(ir);
    end
    
    %     The final point
    
    for in= 1:n_mnr
        line (in,l_size+1)= R (end_coord(1),end_coord(2),end_coord(3),in);
    end
    x_plot(l_size+1) = l_size+1;
    
    % Rescaling the line
    
    for in = 1:n_mnr
        for ip = 1:length(x_plot)
            line_new(in,ip) = weight(in)*(line(in,ip)-polmin(in))/l_length(in);
            
        end
    end
    % Finding intersections
    clear x_inter inter_point intervalue
    
    n_intervalue  = 0;
    for in = 1:n_mnr
        n_intervalue  = n_intervalue  + (in-1);
    end
    
    x_inter_store = cell(1,n_intervalue);
    loop_round = 1;
    k = 0;
    involved = zeros(n_mnr,n_mnr-1);
    while n_mnr-loop_round > 0
        
        for j = 1:n_mnr-loop_round
            k = k+1;
            involved(loop_round,j)= k;  %Needed for isovalue calculation
            involved(loop_round+j,n_mnr-j)= k;
            if ~isempty(find(diff(sign(line_new(loop_round+j,:)-line_new(loop_round,:))), 1))
                x_inter_store{k} = find(diff(sign(line_new(loop_round+j,:)-line_new(loop_round,:))));
                
                
                for i = 1:length(cell2mat(x_inter_store(k)))
                    x_inter = cell2mat(x_inter_store(k));
                    p(1,1)= line_new(loop_round,x_inter(i));
                    p(1,2)= line_new(loop_round,1+x_inter(i));
                    p(2,1)= line_new(loop_round+j,x_inter(i));
                    p(2,2)= line_new(loop_round+j,1+ x_inter(i));
                    inter_point(k,i) = (p(1,1)*(p(2,2)-p(2,1))-p(2,1)*(p(1,2)-p(1,1)))/(p(2,2)-p(2,1)-p(1,2)+p(1,1));
                end
                ipts = inter_point(k,:);
                if ~isempty(ipts(ipts > 0.05 & ipts< 0.95))
                    intervalue(k) = max(ipts(ipts < 0.95));
                else
                    intervalue(k) = max(inter_point(k,:));
                end
            end
        end
        loop_round = loop_round +1;
    end
    
    in_mat = zeros(n_mnr,n_mnr-1);
    
    for k = 1:n_mnr-1
        for in = 1:n_mnr
            in_mat(in,k)= intervalue(involved(in,k));
            isovalue_s(in) = max(in_mat(in,:));
            
        end
    end
    
    %Closing any gaps
    if n_mnr ==3 && range(isovalue_s) > 1e-4
        for in = 1:n_mnr
            r_row = in_mat(in,:);
            no_2(in) = max(r_row(r_row <max(r_row)));
        end
        
        for in = 1:n_mnr-1
            if isovalue_s(in)== no_2(in+1)
                start_mat = find(diff(sign(isovalue_s(in+1)-line_new(in+1,:))));
                intersect_val = (find(diff(sign(line_new(in,:)-line_new(in+1,:))), 1));
                start_ind = start_mat(abs(start_mat - intersect_val)==min(abs(start_mat - intersect_val)));
                py1= line_new(in+1,start_ind);
                py2= line_new(in+1,1+start_ind);
                px1= start_ind;
                px2= start_ind+1;
                y_int = isovalue_s(in+1);
                x_int = px1 + ((y_int-py1)*(px2-px1))/(py2-py1);
                
                
                py1= line_new(in,start_ind);
                py2= line_new(in,1+start_ind);
                px1= start_ind;
                px2= start_ind+1;
                y_new = py1 + ((x_int-px1)*(py2-py1))/(px2-px1);
                isovalue_s(in) = y_new;
            end
        end
        
        for in = 2:n_mnr
            if isovalue_s(in)== no_2(in-1)
                start_mat = find(diff(sign(isovalue_s(in-1)-line_new(in-1,:))));
                intersect_val = (find(diff(sign(line_new(in,:)-line_new(in-1,:))), 1));
                start_ind = start_mat(abs(start_mat - intersect_val)==min(abs(start_mat - intersect_val)));
                
                py1= line_new(in-1,start_ind);
                py2= line_new(in-1,1+start_ind);
                px1= start_ind;
                
                px2= start_ind+1;
                y_int = isovalue_s(in-1);
                x_int = px1 + ((y_int-py1)*(px2-px1))/(py2-py1);
                
                
                py1= line_new(in,start_ind);
                py2= line_new(in,1+start_ind);
                px1= start_ind;
                px2= start_ind+1;
                y_new = py1 + ((x_int-px1)*(py2-py1))/(px2-px1);
                isovalue_s(in) = y_new;
            end
        end
        
        if isovalue_s(n_mnr) == no_2(1)
            start_mat = find(diff(sign(isovalue_s(1)-line_new(1,:))));
            intersect_val = (find(diff(sign(line_new(n_mnr,:)-line_new(1,:))), 1));
            start_ind = start_mat(abs(start_mat - intersect_val)==min(abs(start_mat - intersect_val)));
            
            py1= line_new(1,start_ind);
            py2= line_new(1,1+start_ind);
            px1= start_ind;
            
            px2= start_ind+1;
            y_int = isovalue_s(1);
            x_int = px1 + ((y_int-py1)*(px2-px1))/(py2-py1);
            
            
            py1= line_new(n_mnr,start_ind);
            py2= line_new(n_mnr,1+start_ind);
            px1= start_ind;
            px2= start_ind+1;
            y_new = py1 + ((x_int-px1)*(py2-py1))/(px2-px1);
            isovalue_s(n_mnr) = y_new;
        end
        
        if isovalue_s(1) == no_2(n_mnr)
            start_mat = find(diff(sign(isovalue_s(n_mnr)-line_new(n_mnr,:))));
            intersect_val = (find(diff(sign(line_new(1,:)-line_new(n_mnr,:))), 1));
            start_ind = start_mat(abs(start_mat - intersect_val)==min(abs(start_mat - intersect_val)));
            
            py1= line_new(n_mnr,start_ind);
            py2= line_new(n_mnr,1+start_ind);
            px1= start_ind;
            px2= start_ind+1;
            y_int = isovalue_s(n_mnr);
            x_int = px1 + ((y_int-py1)*(px2-px1))/(py2-py1);
            
            
            py1= line_new(1,start_ind);
            py2= line_new(1,1+start_ind);
            px1= start_ind;
            px2= start_ind+1;
            y_new = py1 + ((x_int-px1)*(py2-py1))/(px2-px1);
            isovalue_s(1) = y_new;
        end
    end
    % Convert the isovalue
    
    for in=1:n_mnr
        isovalue(in) = (isovalue_s(in)*l_length(in))/weight(in) +polmin(in);
    end
end

%% 3-D Visualisation


% Drawing maps
polmaxf = zeros(1,n_mnr);
cn = zeros (1,ncolour);

for j = 1:ncolour
    cn(j)=(10^n_dp);
end

for in = 1:n_mnr
    face1(:,in) = reshape(squeeze(R(1,:,:,in)),[],1);
    face2(:,in) = reshape(squeeze(R(:,1,:,in)),[],1);
    face3(:,in) = reshape(squeeze(R(:,:,1,in)),[],1);
    face4(:,in) = reshape(squeeze(R(grid(1)+1,:,:,in)),[],1);
    face5(:,in) = reshape(squeeze(R(:,grid(2)+1,:,in)),[],1);
    face6(:,in) = reshape(squeeze(R(:,:,grid(3)+1,in)),[],1);
    face_data(:,in) = [face1(:,in); face2(:,in); face3(:,in); face4(:,in); face5(:,in); face6(:,in)];
    polmaxf(in) = max(face_data(:,in)); %Max Density value on face for polymer i
    cn(map_choice(in)) = 1+ ceil((10^n_dp)*polmaxf(in)-((10^n_dp)*isovalue(in))); %Effective colourmap range (+5 to buffer)
    newisovalue(in) = in + isovalue(in) - 1;
    mono_label(in) = char(in+'A'-1);
    titles(in) = {[mono_label(in) ' block Density Profile']};
end


% low = low fraction, i.e. light

colour_low = zeros (ncolour,3);
colour_low(1,:) = [0,0.7,0.9]; %blue
colour_low(2,:) = [0.9,0,0]; %red
colour_low(3,:) = [0,0.9,0.2]; %green
colour_low(4,:) = [1,1,0]; %yellow
colour_low(5,:) = [0.5,0,1]; %purple
colour_low(6,:) = [1,0,1]; %pink
colour_low(7,:) = [0.75,0.75,0.75]; %grey
colour_low(8,:) = [1,1,1];

% high = high fraction, i.e. dark

colour_high = zeros (ncolour,3);
colour_high(1,:) = [0,0,0.4]; %blue
colour_high(2,:) = [0.4,0,0]; %red
colour_high(3,:) = [0,0.4,0]; %green
colour_high(4,:) = [0.4,0.4,0]; %yellow
colour_high(5,:) = [0.15,0,0.30]; %purple
colour_high(6,:) = [0.25,0,0.25]; %pink
colour_high(7,:) = [0.3,0.3,0.3]; %grey
colour_high(8,:) = [1,1,1];

colourpad(:,:,1) = colour_low;
colourpad(:,:,2) = colour_high;
outcolor = colourpad(:,:,1);

map_store = cell(1,ncolour);



for in = 1:ncolour
    temp_map = zeros(cn(in),3);
    temp_map(:,1) = linspace(colourpad(in,1,1),colourpad(in,1,2),cn(in)); %Red
    temp_map(:,2) = linspace(colourpad(in,2,1),colourpad(in,2,2),cn(in)); %Green
    temp_map(:,3) = linspace(colourpad(in,3,1),colourpad(in,3,2),cn(in)); %Blue
    map_store{in}=temp_map;
end

for in = mono_disp
    
    figure(in)
    title(titles(in))
    data = R(:,:,:,in);
    p1 = patch(isosurface(x,y,z,data,isovalue(in)), ...
        'FaceColor',outcolor(map_choice(in),:),'EdgeColor','none','FaceAlpha',opacity(in,1));
    p2 = patch(isocaps(x,y,z,data,isovalue(in)), ...
        'FaceColor','interp','EdgeColor','none','FaceAlpha',opacity(in,2));
    colormap(cell2mat(map_store(map_choice(in))))
    
    if isovalue(in) < max(face_data(:,in))
                      
        cblabelstart = isovalue(in);
        cblabelend = polmaxa(in);
        if cblabelend - cblabelstart > 0.1
            cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),2);
        else
            cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),3);
        end
        l_lngth = linspace(cblabelstart,cblabelend,10);
        cbh = colorbar;
        set(cbh,'ylim',[cblabelstart cblabelend],'ytick',l_lngth,'Yticklabel',cblabel(in,:))
        
        title1 = {['\fontsize{13}\Phi' '_' mono_label(in)]};
        title(cbh,title1)
            
                
    else
        cell1 = {'\fontsize{14}\Phi = ' , num2str(round(isovalue(in),2))};
        text_disp = {strjoin(cell1)};
        %[s1,s2,s3] = ind2sub([plot_grid n_mnr],find(abs((R(:,:,:,in) - isovalue(in)))<0.001,1));
        text( x(grid(1)+1,1,round(grid(3)/2)), y(grid(1)+1,1,round(grid(3)/2)), z(grid(1)+1,1,round(grid(3)/2)),[text_disp])
    end
    
    if strcmp(type,'hexagonal') == 1
        for i =1:2
            size = (grid(1)+1)*(grid(2)+1)*(grid(3)+1);
            coord_set = zeros(size,3);
            counter = 0;
            rotangle = 2*pi/3;
            
            for iz = 1:grid(3)+1
                for iy = 1:grid(2)+1
                    for ix = 1:grid(1)+1
                        counter = counter +1;
                        coord_set(counter,1) = x(ix,iy,iz) ;
                        coord_set(counter,2) = y(ix,iy,iz) ;
                        coord_set(counter,3) = z(ix,iy,iz) ;
                    end
                end
            end
            
            coord_set = coord_set*[cos(rotangle),sin(rotangle),0;-sin(rotangle),cos(rotangle),0;0,0,1];
            
            counter = 0;
            for iz = 1:grid(3)+1
                for iy = 1:grid(2)+1
                    for ix = 1:grid(1)+1
                        counter = counter +1;
                        x(ix,iy,iz) = coord_set(counter,1) ;
                        y(ix,iy,iz) = coord_set(counter,2) ;
                        z(ix,iy,iz) = coord_set(counter,3) ;
                    end
                end
            end
            
            figure(in)
            data = R(:,:,:,in);
            p1 = patch(isosurface(x,y,z,data,isovalue(in)), ...
                'FaceColor',outcolor(map_choice(in),:),'EdgeColor','none','FaceAlpha',opacity(in,1));
            p2 = patch(isocaps(x,y,z,data,isovalue(in)), ...
                'FaceColor','interp','EdgeColor','none','FaceAlpha',opacity(in,2));
            
        end
    end
    
    if dim == 3
        view(3);                        %Sets 3-D view
    elseif dim == 2
        view(2);                        %Sets 2-D view
    elseif dim == 1
        view(2);                        %Sets 2-D view
    end
    
    
    axis equal;                     %Equates the aspect ratio for each axis
    axis vis3d;                     %Freezes aspect ratio (allowing rotation)
    %axis tight;                     %Snaps the axis to the data set
    
    draw_lattice(cell_d,angle,thick,box_clr)
    
    set(gcf,'Renderer','zbuffer');
    %     set(p2,'AmbientStrength',.6);
    %     set(p1,'AmbientStrength',.5);
    
    % isonormals(data,p1);
    % lightangle(45,30);
    % lighting phong
    
end

if ~isempty(comp_disp)
    draw_all (newisovalue,opacity,n_mnr,R,x,y,z,n_dp,map_store,outcolor,cn,type,grid,comp_disp,map_choice,colourpad,isovalue,dim,polmaxa)
    draw_lattice(cell_d,angle,thick,box_clr)
    if n_mnr>4
        figure(n_mnr+1)
        cbh=colorbar;
        c = 0;
        for in = 1:n_mnr
            c = c+1;
            leftlabel(c) = round(isovalue(in)+in -1,2);
            leftlabel2(c) = round(isovalue(in),2);
            c = c+1;
            leftlabel(c) = round(polmaxa(in)+in -1,2);
            leftlabel2(c) = round(polmaxa(in),2);
        end
        
        rightlabel = [0:0.1:1 0.1:0.1:1 0.1:0.1:1];
        h1=axes('position',get(cbh,'position'),'color','none','ylim',[isovalue(1)-0.01,polmaxa(n_mnr)+2.01],'ytick',leftlabel,'yticklabel',leftlabel2,'xtick',[]);
        set(cbh,'YTick',0:0.1:3,'Yticklabel',rightlabel)
        title(cbh,'\fontsize{14}\Phi')
    else
        if n_mnr < 4
        cb_pos(1,:)= [.68 .11 .04 .815]; %x,y,width,length
        cb_pos(2,:)= [.79 .11 .04 .815]; %x,y,width,length
        cb_pos(3,:)= [.90 .11 .04 .815]; %x,y,width,length
        else
        cb_pos(1,:)= [.57 .11 .04 .815]; %x,y,width,length 
        cb_pos(2,:)= [.68 .11 .04 .815]; %x,y,width,length
        cb_pos(3,:)= [.79 .11 .04 .815]; %x,y,width,length
        cb_pos(4,:)= [.90 .11 .04 .815]; %x,y,width,length    
        end
        
        figure(n_mnr+1)
        ax(n_mnr+1) = gca;
        
        clear cblabel
        for in = 1:n_mnr
            c=0;
            
            figure(n_mnr+1)
            
            ax(in) = axes;
            colormap(ax(in),cell2mat(map_store(map_choice(in))))
            ax(in)= gca; 
            ax(in).Visible = 'off';
            ax(in).XTick = [];
            ax(in).YTick = [];
            ax(in).ZTick = [];
            
            if isovalue(in) < max(face_data(:,in))
                
            cblabelstart = isovalue(in);
            cblabelend = polmaxa(in);
            if cblabelend - cblabelstart > 0.1
            cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),2);
            else
            cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),3);    
            end
            l_lngth = linspace(0,1,10);
            cb(in) = colorbar(ax(in),'Position',cb_pos(in,:));
            set(cb(in),'ytick',l_lngth,'Yticklabel',cblabel(in,:))
            %set(cb(in),'ylim',[isovalue(in)-0.01,polmaxa(in)+0.01])
            title1 = {['\fontsize{13}\Phi' '_' mono_label(in)]};
            title(cb(in),title1)
            end
        end
        
        linkprop(ax,{'view'});
        if n_mnr ~= 4
        set(ax,'Position',[.05 .11 .55 .815]); %x,y,width,length 
        else
        set(ax,'Position',[.05 .11 .44 .815]); %x,y,width,length  
        end
        
        if dim == 3
            view(3);                        %Sets the view to 3-D
        else
            view(2);                        %Sets the view to 2-D
        end
        % ax(n_mnr+1) = gca;
        % axis vis3d;                     %Freezes aspect ratio (allowing rotation)
    end
end

%% Plotting the lines

if linedraw == 1
    figure(n_mnr+ 2)
    
    linestyles(1) = {'-.'};
    linestyles(2) = {':'};
    linestyles(3) = {'--'};
    for in = 1:n_mnr
        plot (x_plot,line_new(in,:),'color',outcolor(map_choice(in),:))
        hold on
    end
    
    clear yy
    for in = 1:n_mnr
        yy(in,:) = ones(1,length(x_plot))*isovalue_s(in);
    end
    
    for in = 1:n_mnr
        figure (n_mnr+2)
        hold on
        plot(x_plot,yy(in,:),'color',outcolor(map_choice(in),:),'LineWidth',1+n_mnr-in,'LineStyle',cell2mat(linestyles(mod(in,3)+1)))
    end
    
    xlabel('Line Index')
    ylabel('Rescaled Volume Fraction')
    
    for in=1:n_mnr
        mono_label(in) = char(in+'A'-1);
        legend_labels1(in) = {[mono_label(in) ' block Rescaled Density']};
        legend_labels2(in) = {[mono_label(in) ' block Rescaled Isovalue']};
    end
    
    
    legend([legend_labels1 legend_labels2])
    legend('Location','best')
    
end

if ~isempty(drawscatter)
    % Scattering
    
    % Pre-setting the matrix
    
    count = 0;
    
    for i = 1:3
        b(i) = (2*pi)/cell_d(i); %basis vectors
    end
    
    
    F_store = zeros(length(h_set)*length(k_set)*length(l_set),5);
    
    for h = h_set
        for k = k_set
            for l = l_set
                count = count + 1;
                F_store(count,1) = h;
                F_store(count,2) = k;
                F_store(count,3) = l;
            end
        end
    end
    
    % Calculating the structure factor
    x_index_label = cell(length(F_store),1);
    x_label = cell(length(F_store),1);
    
    I = zeros(1,length(F_store));
    q = zeros(1,length(F_store));
    x_index = zeros(1,length(F_store));
    
    for i_f = 1:length(F_store)
        h = F_store(i_f,1);
        k = F_store(i_f,2);
        l = F_store(i_f,3);
        F_sum = 0;
        for in = drawscatter
            x_s = zeros(grid);
            y_s = zeros(grid);
            z_s = zeros(grid);
            
            for iz=1:grid(3)
                for iy=1:grid(2)
                    for ix=1:grid(1)
                        x_s(ix,iy,iz) = x(ix,iy,iz)/cell_d(1);
                        y_s(ix,iy,iz) = y(ix,iy,iz)/cell_d(2);
                        z_s(ix,iy,iz) = z(ix,iy,iz)/cell_d(3);
                        F_sum = F_sum + R(ix,iy,iz,in)*exp(2*1i*pi*((h*x_s(ix,iy,iz))+(k*y_s(ix,iy,iz))+(l*z_s(ix,iy,iz))));
                        
                    end
                end
            end
            
        end
        F_store(i_f,4) = F_sum * cell_d(1)*cell_d(2)*cell_d(3);
        F_store(i_f,5) = abs(F_sum)^2;
        I(i_f) = F_store(i_f,5);
        x_index(i_f) = i_f;
        x_index_label(i_f) = {[ h k l ]};
        x_label{i_f}=mat2str(cell2mat(x_index_label(i_f)));
        % x_index_mag(i_f) = h^2 + k^2 + l^2;
        % x_index_3(i_f) = h + k + l;
        q(i_f) = sqrt((b(1)*h)^2 + (b(2)*k)^2 + (b(3)*l)^2);
    end
    
    % Plotting scatterplots
    
    %I(q)
    
    plotmat = sortrows([q;I]');
    
    [q1,~,q_ind] = uniquetol(plotmat(:,1));
    
    plotmat_avg = [q1,accumarray(q_ind,plotmat(:,2),[],@mean)];
    
    q_sort = plotmat_avg(:,1);
    I_sort = plotmat_avg(:,2);
    
    % Making a delta function
    
    t_fcr = 100;
    q2 = linspace(0,max(q_sort),length(q_sort)*t_fcr);
    q_plot = sort([q2 q_sort']);
    I_plot = zeros(length(q_plot),1)+min(I);
    
    step = 1;
    for i = 1:length(q_plot)
        if q_plot(i)== q_sort(step)
            
            I_plot(i) = I_sort(step);
            step = step + 1;
            if step > length(q_sort)
                break
            end
        end
    end
    
    figure(n_mnr+3)
    semilogy(q_plot,I_plot,'.-')
    %plot(q_plot,I_plot,'.-')
    xlabel(['q [' char(197) '^-^1]'])
    ylabel('Intensity')
    
    % vs Miller Indices
    
    t_fcr = 100;
    x_index_plot = linspace(1,length(x_index)*t_fcr,length(x_index)*t_fcr);
    I_plot2 = zeros(length(x_index_plot),1)+min(I);
    for ii = 1:length(x_index)
        I_plot2((ii*t_fcr)-(t_fcr-1)) = I(ii);
    end
    
    figure(n_mnr+4)
    semilogy(x_index_plot,I_plot2,'.-')
    %plot(x_index_plot,I_plot2,'.-')
    ylabel('Intensity')
    xlabel('Miller Indices')
    
    x_label{length(x_label)+1}=' ';
    
    ax = gca;
    %ax.YLim = [0 max(I)];
    ax.XTick = linspace(1,length(x_index)*t_fcr,length(x_index)+1);
    ax.XTickLabel = x_label;
    ax.XTickLabelRotation = 45;
    
end

if ~isempty(inputvec)
    
    figure(n_mnr+5)
    l_size = 0;
    userinput = inputvec;   
    start_coord = [1 1 1];
    end_coord = userinput .* grid + [1 1 1];
    
    dir_vec = end_coord-start_coord;
    step_length = max(abs(dir_vec));
    clear ix iy iz x_plot
    
    for il = 1:step_length
        ix(il) = start_coord(1)+ round((il-1)*(dir_vec(1)/step_length));
        iy(il) = start_coord(2)+ round((il-1)*(dir_vec(2)/step_length));
        iz(il) = start_coord(3)+ round((il-1)*(dir_vec(3)/step_length));
        x_plot(il) = (il-1)/(step_length-1);
        for in= 1:n_mnr
            line_plot (in,il) = R (ix(il),iy(il),iz(il),in);
            
        end
    end
    
    linemarkers(1) = {'s'};
    linemarkers(2) = {'^'};
    linemarkers(3) = {'o'};
    linemarkers(4) = {'p'};
    linemarkers(5) = {'+'};
    linemarkers(6) = {'*'};
    linemarkers(7) = {'d'};
    
    for in = 1:n_mnr
        plot (x_plot,line_plot(in,:),'color',outcolor(map_choice(in),:),'marker',cell2mat(linemarkers(mod(in,3)+1)),'markerfacecolor',outcolor(map_choice(in),:),'markersize',5)
        hold on
        legend_labels3(in) = {[mono_label(in) ' block']};
        
    end
    
    
    title1 = {['Density Variation Along ' '[' num2str(userinput) ']']};
    title(strjoin(title1))
    xlabel('\fontsize{14}r/r\fontsize{14}_m_a_x')
    ylabel('\fontsize{14}\Phi\fontsize{13}(r)')
    %ylim([0 1])
    
    legend([legend_labels3])
    legend('Location','northeastoutside')
    
    
    
end



%out = isovalue
toc

%%
%keyboard
end

% Subfunctions


function draw_all (newisovalue,opacity,n_mnr,R,x,y,z,n_dp,map_store,outcolor,cn,type,grid,comp_disp,map_choice,colourpad,isovalue,dim,polmaxa)

comp_disp = sort(comp_disp);

for in = 1:n_mnr
    face1 = reshape(squeeze(R(1,:,:,in)),[],1);
    face2 = reshape(squeeze(R(:,1,:,in)),[],1);
    face3 = reshape(squeeze(R(:,:,1,in)),[],1);
    face4 = reshape(squeeze(R(grid(1)+1,:,:,in)),[],1);
    face5 = reshape(squeeze(R(:,grid(2)+1,:,in)),[],1);
    face6 = reshape(squeeze(R(:,:,grid(3)+1,in)),[],1);
    face_data(:,in) = [face1; face2; face3; face4; face5; face6];
end


fill = zeros(1,n_mnr-1);

c = 0;
for in = comp_disp(1:end-1)
    c = c +1;
    
    if in~= n_mnr
        fill(in) = (newisovalue(comp_disp(c+1))-newisovalue(in))*(10^n_dp) - cn(map_choice(in));
    end
    
end

temp_map = ones(2*(10^n_dp),3);

fillmap_store = cell(1,n_mnr-1);

c = 0;

for in = comp_disp(1:end-1)
    c = c +1;
    temp_fillmap = zeros(round(fill(in)),3);
    temp_fillmap(:,1) = linspace(colourpad(map_choice(comp_disp(c)),1,2),colourpad(map_choice(comp_disp(c+1)),1,1),round(fill(in))); %Red
    temp_fillmap(:,2) = linspace(colourpad(map_choice(comp_disp(c)),2,2),colourpad(map_choice(comp_disp(c+1)),2,1),round(fill(in))); %Green
    temp_fillmap(:,3) = linspace(colourpad(map_choice(comp_disp(c)),3,2),colourpad(map_choice(comp_disp(c+1)),3,1),round(fill(in))); %Blue
    fillmap_store{in} = temp_fillmap;
end

% for in = comp_disp(1:end-1)
%    temp_fillmap = zeros(round(fill(in)),3);
%     for ifl =1:fill(in)
%         temp_fillmap(ifl,:) = temp_map(ifl,:);
%     end
%     fillmap_store{in} = temp_fillmap;
% end


newmap = [];
for in = comp_disp
    if isovalue(in) < max(face_data(:,in))
        if in == n_mnr
            newmap = [newmap;map_store(map_choice(in))];
        else
            newmap = [newmap;map_store(map_choice(in));fillmap_store(in)];
        end
    end
end

% Simultaneous Visualisation
D = R;
figure(n_mnr+1)
title('Composite Density Profile')
for in = comp_disp
    D(:,:,:,in) = R(:,:,:,in) +in -1;
    data = (D(:,:,:,in));
    p1 = patch(isosurface(x,y,z,data,newisovalue(in)), ...
        'FaceColor',outcolor(map_choice(in),:),'EdgeColor','none','FaceAlpha',opacity(in,1));
    p2 = patch(isocaps(x,y,z,data,newisovalue(in)), ...
        'FaceColor','interp','EdgeColor','none','FaceAlpha',opacity(in,2));
    
    if dim == 3
        view(3);                        %Sets the view to 3-D
    else
        view(2);                        %Sets the view to 2-D
    end
    axis equal;                     %Equates the aspect ratio for each axis
    axis vis3d;                     %Freezes aspect ratio (allowing rotation)
    %axis tight;                     %Snaps the axis to the data set
    
    if strcmp(type,'hexagonal') == 1
        size = (grid(1)+1)*(grid(2)+1)*(grid(3)+1);
        coord_set = zeros(size,3);
        
        for i =1:2
            
            counter = 0;
            rotangle = 2*pi/3;
            
            for iz = 1:grid(3)+1
                for iy = 1:grid(2)+1
                    for ix = 1:grid(1)+1
                        counter = counter +1;
                        coord_set(counter,1) = x(ix,iy,iz) ;
                        coord_set(counter,2) = y(ix,iy,iz) ;
                        coord_set(counter,3) = z(ix,iy,iz) ;
                    end
                end
            end
            
            coord_set = coord_set*[cos(rotangle),sin(rotangle),0;-sin(rotangle),cos(rotangle),0;0,0,1];
            
            counter = 0;
            for iz = 1:grid(3)+1
                for iy = 1:grid(2)+1
                    for ix = 1:grid(1)+1
                        counter = counter +1;
                        x(ix,iy,iz) = coord_set(counter,1) ;
                        y(ix,iy,iz) = coord_set(counter,2) ;
                        z(ix,iy,iz) = coord_set(counter,3) ;
                    end
                end
            end
            
            data = D(:,:,:,in);
            p1 = patch(isosurface(x,y,z,data,newisovalue(in)), ...
                'FaceColor',outcolor(map_choice(in),:),'EdgeColor','none','FaceAlpha',opacity(in,1));
            p2 = patch(isocaps(x,y,z,data,newisovalue(in)), ...
                'FaceColor','interp','EdgeColor','none','FaceAlpha',opacity(in,2));
            set(gcf,'Renderer','zbuffer')
        end
    end
    
    
end


colormap(cell2mat(newmap))

% cbh= colorbar ;
% set(cbh,'YTick',0:0.1:3)
%
% title(cbh,'\fontsize{14}\Phi')
% cbh.Location = 'eastoutside';
% cbh.Color = [1 0 1];
end

function draw_lattice(cell_d,angle,thick,box_clr)
%%
axes = gca;
x_pos(1) = 0;
x_pos(2) = cos(angle(3))*cell_d(2);
x_pos(3) = cell_d(1)+ cos(angle(3))*cell_d(2);
x_pos(4) = cell_d(1);
x_pos(5) = 0 + cos(angle(1))*cell_d(3);
x_pos(6) = cos(angle(3))*cell_d(2) + cos(angle(1))*cell_d(3);
x_pos(7) = cos(angle(3))*cell_d(2) + cos(angle(1))*cell_d(3) + cell_d(1);
x_pos(8) = cell_d(1) + cos(angle(1))*cell_d(3);

x_start = x_pos;
for i =9:12
    x_start(i)=x_pos(i-8);
end
x_end = zeros(length(x_start),1)';

for i =1:4
    if i == 4
        x_end(i)=x_start(1);
    else
        x_end(i)=x_start(i+1);
    end
end


for i =5:8
    if i == 8
        x_end(i)=x_start(5);
    else
        x_end(i)=x_start(i+1);
    end
end

for i =9:12
    x_end(i)=x_start(i-4);
end

X1 = [x_start;x_end];

y_pos(1) = 0;
y_pos(2) = sin(angle(3))*cell_d(2);
y_pos(3) = sin(angle(3))*cell_d(2);
y_pos(4) = 0;
y_pos(5) = 0 + cos(angle(2))*cell_d(3);
y_pos(6) = sin(angle(3))*cell_d(2) + cos(angle(2))*cell_d(3);
y_pos(7) = sin(angle(3))*cell_d(2) + cos(angle(2))*cell_d(3);
y_pos(8) = 0 + cos(angle(2))*cell_d(3);

y_start = y_pos;
for i =9:12
    y_start(i)=y_pos(i-8);
end
y_end = zeros(length(y_start),1)';

for i =1:4
    if i == 4
        y_end(i)=y_start(1);
    else
        y_end(i)=y_start(i+1);
    end
end

for i =5:8
    if i == 8
        y_end(i)=y_start(5);
    else
        y_end(i)=y_start(i+1);
    end
end

for i =9:12
    y_end(i)=y_start(i-4);
end

Y1 = [y_start;y_end];

z_pos(1) = 0;
z_pos(2) = 0;
z_pos(3) = 0;
z_pos(4) = 0;
z_pos(5) = cell_d(3)*sin(angle(1))*sin(angle(2));
z_pos(6) = cell_d(3)*sin(angle(1))*sin(angle(2));
z_pos(7) = cell_d(3)*sin(angle(1))*sin(angle(2));
z_pos(8) = cell_d(3)*sin(angle(1))*sin(angle(2));

z_start = z_pos;
for i =9:12
    z_start(i)=z_pos(i-8);
end
z_end = zeros(length(z_start),1)';

for i =1:4
    if i == 4
        z_end(i)=z_start(1);
    else
        z_end(i)=z_start(i+1);
    end
end


for i =5:8
    if i == 8
        z_end(i)=z_start(5);
    else
        z_end(i)=z_start(i+1);
    end
end

for i =9:12
    z_end(i)=z_start(i-4);
end

Z1 = [z_start;z_end];

line(X1,Y1,Z1,'color',box_clr,'LineStyle','-','LineWidth',thick)

end