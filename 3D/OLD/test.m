

[x,y,z,v] = flow;

levellist = linspace(-10,2,7);

for i = 1:length(levellist)
    level = levellist(i);
    p = patch(isosurface(x,y,z,v,level));
    p.FaceVertexCData = level;
    p.FaceColor = 'flat';
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.3;
end
view(3)