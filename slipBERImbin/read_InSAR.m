function [ locs_InSAR, d_InSAR, los_vector_InSAR ] = read_InSAR( input_filename )

% Ruth Amey March 2018 (making previous, silly codes less silly)

fin = fopen(deblank(input_filename),'r');
displ_input = fscanf(fin,'%lf',[6 inf]);
displ_input = displ_input';
fclose(fin);
disp(['  processing file: ' input_filename]) ;

locs_InSAR(1,:)=displ_input(:,1) ;
locs_InSAR(2,:)=displ_input(:,2) ;
locs_InSAR(3,:)=0 ;

d_InSAR(:)=displ_input(:,3) ;

los_vector_InSAR(1,:)=displ_input(:,4) ;
los_vector_InSAR(2,:)=displ_input(:,5) ;
los_vector_InSAR(3,:)=displ_input(:,6) ;

end

