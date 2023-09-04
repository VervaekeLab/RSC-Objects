function extractedFrames = getSelectedFrames(fileName,selected_frames)


% Open file for reading at the beginning
fid = fopen(fileName,'r','b');

% Get header info
headerInfo = sb.database.load.getHeaderInfo(fid);

% Use the Little Endian machine format ordering for reading bytes
endianType = 'ieee-le';

%% Reading images
% Following the header, each images is stored and aligned on a 8192 bytes
% boundary (when no metadata is included)
imageOffset = 8192;

bitstr = '';

switch headerInfo.ImageBitDepthReal
    case 8
        bitstr = 'uint8';
    case {10,12,14,16}
        bitstr = 'uint16';
end
if isempty(bitstr)
    error('Unsupported bit depth');
end

% Initialise variable for the number of frames read
nread = 1;
        
numPixels = headerInfo.ImageWidth * headerInfo.ImageHeight;

% Start the read loop
count = 1;

% Allocate variable 
extractedFrames = zeros(headerInfo.ImageHeight,headerInfo.ImageWidth,length(selected_frames));

for f = 1:length(selected_frames)
   
    % Find image current position
    nread = selected_frames(f);
    
    % Move to correct position in file
    fseek(fid, imageOffset + nread * headerInfo.TrueImageSize, 'bof');
    
    % Read, interpret and convert the data dependent on its format
    MonoColVec = fread(fid,numPixels,bitstr,endianType);
    
    img = reshape(MonoColVec,headerInfo.ImageWidth,headerInfo.ImageHeight)';
    extractedFrames(:,:,f) = img;
    
    
end

fclose(fid);

end