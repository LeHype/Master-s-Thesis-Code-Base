function [Position] = GetPosition_fromCrane()
cr = Crane3D;
Position = get(cr,'Position');
end