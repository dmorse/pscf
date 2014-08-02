% 2 monomer, 3 dimension
[v1,v2] = textread('di/di-o70_net2.dat','%f%f','delimiter','\n');
if v1(1) == 1
    cell(1:3) = v1(2);   % 1 parameter system, e.g. Gyroid
elseif v1(1) == 3
    cell(1:3) = v1(2:4); % 3 parameter system, e.g. O70
end
counter = 1 + v1(1);
grid = v1(counter+1 : counter+3);
counter = counter + 3;
clear x y z a b c b2;
for iz=1:grid(3),
    for iy=1:grid(2),
        for ix=1:grid(1),
            counter = counter + 1;
            x(ix,iy,iz) = cell(1) * (ix-1)/grid(1);
            y(ix,iy,iz) = cell(2) * (iy-1)/grid(2);
            z(ix,iy,iz) = cell(3) * (iz-1)/grid(3);
            a(ix,iy,iz) = v1(counter);
            b(ix,iy,iz) = v2(counter);
        end
    end
end
    
for iz=1:grid(3),
    for iy=1:grid(2),
        for ix=1:grid(1),
            if (iz > grid(3)/2)
                b2(ix,iy,iz) = b(ix,iy,iz-grid(3)/2);
            else
                b2(ix,iy,iz) = b(ix,iy,iz+grid(3)/2);
            end
        end
    end
end

box on;
data = smooth3(a,'box',5);
isovalue =.5;
p1 = patch(isosurface(x,y,z,data,isovalue), ...
     'FaceColor','blue','EdgeColor','none');
p2 = patch(isocaps(x,y,z,data,isovalue), ...
     'FaceColor','interp','EdgeColor','none');
isonormals(data,p1);
axis equal;
view(3); axis vis3d tight
lightangle(45,30);
set(gcf,'Renderer','zbuffer'); lighting phong
set(p2,'AmbientStrength',.6); 
set(p1,'AmbientStrength',.5); 
%set(p1,'SpecularColorReflectance',1,'SpecularExponent',50)