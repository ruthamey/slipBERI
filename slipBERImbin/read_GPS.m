function [ locs_GPS, d_GPS, los_vector_GPS, sigma_gps ] = read_GPS( input_filename )

% Ruth Amey March 2018 (making previous, silly codes less silly)

fin = fopen(deblank(input_filename),'r');
displ_input = fscanf(fin,'%lf',[7 inf]);
displ_input = displ_input';
fclose(fin);
disp(['  processing file: ' input_filename]) ;

locs_GPS(1,:)=displ_input(:,1);
locs_GPS(2,:)=displ_input(:,2);
locs_GPS(3,:)=0 ;

d_GPS(:)=displ_input(:,3) ;

los_vector_GPS(1,:)=displ_input(:,4) ;
los_vector_GPS(2,:)=displ_input(:,5) ;
los_vector_GPS(3,:)=displ_input(:,6) ;

sigma_gps = displ_input(:,7);

end

