function Kmeans_results=LEiDA_EEG_kmeans(data, thisK, gpu)

    tic
    msg=sprintf('Starting kmeans for k=%i.', thisK);
    disp(msg)
    if gpu==1
        data=gpuArray(data);
    end
    
    [IDX, C, SUMD, D]=kmeans(data, thisK, 'distance', 'sqeuclidean', 'Replicates', 50, 'MaxIter', 200, 'Display', 'off');
    [~, ind_sort]=sort(hist(IDX,1:thisK),'descend');
    [~,idx_sort]=sort(ind_sort,'ascend');
    Kmeans_results.IDX=idx_sort(IDX);   % Cluster time course - numeric collumn vectos
    Kmeans_results.C=C(ind_sort,:);       % Cluster centroids (FC patterns)
    Kmeans_results.SUMD=SUMD(ind_sort); % Within-cluster sums of point-to-centroid distances
    Kmeans_results.D=D(ind_sort);       % Distance from each point to every centroid
    clear data
    elapsedTime=toc;
    msg=sprintf('Finished kmeans for k=%i in %d seconds.', thisK, round(elapsedTime));
    disp(msg)
end