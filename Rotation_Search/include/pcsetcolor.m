function pcsetcolor( cloud, color )

if color=='r'
    cloud.Color = uint8(repmat([255 0 0],cloud.Count,1));
elseif color=='m'
    cloud.Color = uint8(repmat([255 0 255],cloud.Count,1));

elseif color=='b'
    cloud.Color = uint8(repmat([0 0 255],cloud.Count,1));    
elseif color=='grey'
    cloud.Color = uint8(repmat([180 180 180],cloud.Count,1));    
end

end

