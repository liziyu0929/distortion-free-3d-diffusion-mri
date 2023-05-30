function alpha = NRMSE(I_ref,I,w)

I = abs(I);
I_ref = abs(I_ref);



if nargin > 2
I = I.*w;
I_ref = I_ref.*w;
end

I = I/norm(I(:),'fro');
I_ref = I_ref/norm(I_ref(:),'fro');

% I = I/max(I(:));
% I_ref = I_ref/max(I_ref(:));


diff = I - I_ref;

sosdiff = sqrt(diff(:)'*diff(:));
sosref = sqrt(I_ref(:)'*I_ref(:));
alpha = sosdiff/sosref;

end