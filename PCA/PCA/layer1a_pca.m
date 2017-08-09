%% set parameters
psize = 8;         % size of an image patch (psize*psize)
D0 = psize*psize;   % dimensionality of layer-1 PCA inputs
d1 = D0 - 1;        % dimensionality of layer-1 ICA inputs

%% initialize
load IMAGES;        % load images
M = length(IMAGES); % number of images
COV = zeros(D0,D0); % the covariance matrix

%% main loop
N = 5000;               % number of image patches to load in one batch
patches = zeros(D0,N);  % image patches loaded in one batch
counter = 0;            % to count the total number of image patches loaded
pointer = 0;            % to count the number of image patches loaded in the current batch
for m = 1:10
    % load one image
    fprintf(1, 'processing image %d ... \n', m);
    this_image = IMAGESr(:,:,m); 
    %this_image=reshape(fea(m*160,:),64,64);%此处为曾经尝试利用人脸图像。
    [width, height] = size(this_image);
   % colormap(gray);
   % subplot(1,5,m);%显示图像
   % imagesc(this_image);
   % axis off;

    % extract image patches from this image %提取图像子块
    for r = 1:width-psize+1
        for c = 1:height-psize+1
            pointer = pointer+1;
            patches(:,pointer) = reshape(this_image(r:r+psize-1,c:c+psize-1),D0,1);
            if (pointer == N)
                pointer = 0;
                counter = counter + N;
                patches = patches - repmat(mean(patches),D0,1);
                COV = COV + patches*patches';
            end;
        end;
    end;
 end;
% if there are image patches left un-processed, add their covariance matrix
if (pointer > 0)
    patches = patches(:,1:pointer);
    counter = counter + pointer;
    patches = patches - repmat(mean(patches),D0,1);
    COV = COV + patches*patches';
end;

%% calculate the dewhiteningMatrix and whiteningMatrix  %计算whitenmatrix
C = COV/(counter-1);
[E, L] = eig(C);
[Vals,order] = sort(diag(L),'descend');
E = E(:,order(1:d1));
L = diag(Vals(1:d1).^(-0.5));%此处为PCA算法的一种方法
whiteningMatrix = L*E';
dewhiteningMatrix = E*L^(-1);
save(['P1-' num2str(D0) '-' num2str(d1) '.mat'], 'dewhiteningMatrix', 'whiteningMatrix', 'Vals');