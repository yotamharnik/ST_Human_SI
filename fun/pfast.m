function p = pfast(p)
% Fisher's (1925) method for combination of independent p-values
% Code adapted from Bailey and Gribskov (1998)

    product=prod(p);
    n=length(p);
    if n<=0
        error('pfast was passed an empty array of p-values')
    elseif n==1
        p = product;
        return
    elseif product == 0
        p = 0;
        return
    else
        x = -log(product);
        t=product;
        p=product;
        for i = 1:n-1
            t = t * x / i;
            p = p + t;
        end
    end  