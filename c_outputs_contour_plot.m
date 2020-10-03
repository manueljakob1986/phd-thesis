function c_outputs_contour_plot(matrix, parking_pricing)

% Contour Plot N_s, N_p and parking pricing (linear)
[x,y]=meshgrid(linspace(0,max(matrix(:,2)),500), linspace(0,max(matrix(:,3)),500));
vq = griddata(matrix(2:size(matrix(:,2),1),2),matrix(2:size(matrix(:,3),1),3),parking_pricing,x,y,'linear');
figure
hold on
contourf(x,y,vq)
% title('Contour Plot')
xlabel('Number of searching vehicles')
ylabel('Number of parking vehicles')
c = colorbar;
c.Label.String = 'Parking price (in CHF)';
hold off

% Contour Plot N_s, N_p and parking pricing (natural)
[x,y]=meshgrid(linspace(0,max(matrix(:,2)),500), linspace(0,max(matrix(:,3)),500));
vq = griddata(matrix(2:size(matrix(:,2),1),2),matrix(2:size(matrix(:,3),1),3),parking_pricing,x,y,'natural');
figure
hold on
contourf(x,y,vq)
% title('Contour Plot')
xlabel('Number of searching vehicles')
ylabel('Number of parking vehicles')
c = colorbar;
c.Label.String = 'Parking price (in CHF)';
hold off

% Contour Plot N_s, A and parking pricing (linear)
[x,y]=meshgrid(linspace(0,max(matrix(:,2)),500), linspace(0,max(matrix(:,4)),500));
vq = griddata(matrix(2:size(matrix(:,2),1),2),matrix(2:size(matrix(:,4),1),4),parking_pricing,x,y,'linear');
figure
hold on
contourf(x,y,vq)
% title('Contour Plot')
xlabel('Number of searching vehicles')
ylabel('Number of available parking spaces')
c = colorbar;
c.Label.String = 'Parking price (in CHF)';
hold off

% Contour Plot N_s, A and parking pricing (natural)
[x,y]=meshgrid(linspace(0,max(matrix(:,2)),500), linspace(0,max(matrix(:,4)),500));
vq = griddata(matrix(2:size(matrix(:,2),1),2),matrix(2:size(matrix(:,4),1),4),parking_pricing,x,y,'natural');
figure
hold on
contourf(x,y,vq)
% title('Contour Plot')
xlabel('Number of searching vehicles')
ylabel('Number of available parking spaces')
c = colorbar;
c.Label.String = 'Parking price (in CHF)';
hold off

end




