#' @keywords internal
scale_coords = function(spdata)
{
	min.x <- min(spdata[,1])
	min.y <- min(spdata[,2])
	
	coords.new <- spdata[,1:2]
	coords.new[,1] <- (coords.new[,1] - min.x)+1
	coords.new[,2] <- (coords.new[,2] - min.y)+1
	
	spdata.shift <-cbind(coords.new, spdata[,3])
	
	return(spdata.shift)
	
}
#' @keywords internal
cov_lags = function(nrows, ncols)
{
	jlags <- c(0:(ncols-1), -1*(1:(ncols-1)))
	jlags <- sort(jlags)
	klags <- 0:(nrows-1)
	lags <- expand.grid(jlags, klags)
	return(cbind(lags[,1],lags[,2]))
}
#' @keywords internal
lag_dist = function(locs1, locs2)
{
	n1 <- dim(locs1)[1]
	n2 <- dim(locs2)[1]
	index <- expand.grid(1:n2,1:n1)
	index <- cbind(index[,1], index[,2])
	index <- index[,2:1]
	splags <- c()
	for(i in 1:n1)
	{ 
		lags.d <- cbind(locs1[i,1] - locs2[,1], locs1[i,2] - locs2[,2])
		splags <- rbind(splags, lags.d)
	}
	
	return(cbind(index, splags))
}
#' @keywords internal
chat_jk = function(lagvec, spdata, nrows, ncols)
{	
	N <- dim(spdata)[1]
	cj <- lagvec[1]
	ck <- lagvec[2]
	mysum <- 0
	if(cj >= 0)
	{
		for(u in 1:(ncols-cj))
		{
			for(v in 1:(nrows-ck))
			{	
				loc1 <- which(spdata[,1] == u & spdata[,2] == v)
				z1 <- spdata[loc1,3]
				loc2 <- which(spdata[,1] == u+cj & spdata[,2] == v+ck)
				z2 <- spdata[loc2,3]	
				mysum <- mysum + (z1*z2)
			}	
		}
		chat <- mysum/N	
	}
	
	if(cj < 0)
	{
		for(u in (-cj+1):ncols)
		{
			for(v in 1:(nrows-ck))
			{	
				loc1 <- which(spdata[,1] == u & spdata[,2] == v)
				z1 <- spdata[loc1,3]
				loc2 <- which(spdata[,1] == u+cj & spdata[,2] == v+ck)
				z2 <- spdata[loc2,3]	
				mysum <- mysum + (z1*z2)
			}	
		}
		chat <- mysum/N	
	}
	
	return(chat)
}
#' @keywords internal
est_cov = function(lagsmat, spdata, nrows, ncols)
{
	chat <- apply(lagsmat, 1, FUN = chat_jk, spdata, nrows, ncols)
	return(chat)	
}
#' @keywords internal
cov_complete = function(chat, lags, ret.mat = FALSE)
{
	njlags <- length(unique(lags[,1]))
	klags.computed <- length(unique(lags[,2]))
	nklags <- 2*length(unique(lags[,2])) - 1
	Cmat <- matrix(data = NA, ncol = njlags, nrow = nklags)
	jlags <- unique(lags[,1])
	klags <- unique(lags[,2])
	zero <- which(klags == 0)
	tmp <- klags[-zero]
	klags <- sort( c(klags, -1*tmp) )
	klags <- rev(klags)
	mat.cols <- 1:njlags
	mat.rows <- 1:nklags
	
	for(i in 1:length(chat))
	{
		cj <- lags[i,1]
		ck <- lags[i,2]
		cloc <- which(jlags == cj)
		rloc <- which(klags == ck)
		row.coord <- mat.rows[rloc]
		col.coord <- mat.cols[cloc]
		Cmat[row.coord,col.coord] <-  chat[i]
	}
	
	lags2 <-  c()
	chat2 <-  c()
	for(a in 1:njlags)
	{
		for(b in (klags.computed+1):nklags) 
		{
			a.map <- which(mat.cols == a)
			j.lag <- jlags[a.map]
			b.map <- which(mat.rows == b)
			k.lag <- klags[b.map]
			d.map <- which(klags == -k.lag)
			d <- mat.rows[d.map]
			e.map <- which(jlags == -j.lag)
			e <- mat.cols[e.map]
			Cmat[b,a] <- Cmat[d,e]
			chat2 <- c(chat2, Cmat[d,e])
			lags2 <- rbind(lags2, c(j.lag, k.lag))
		}
	}
	
	chat.full <- cbind(rbind(lags, lags2), c(chat, chat2))
	
	if(ret.mat == T)
	{
		list("covmat" = Cmat, "lags" = rbind(lags, lags2))
	}
	else
	{
		return(chat.full)
	}
}
#' @keywords internal
get_Fourier_freqs = function(nrowss, ncolss)
{
	s <- 1:ncolss
	p <- c( -1*floor((s-1)/2 ), floor(s/2) )
	p <- sort( unique(p) )
	om1 <- 2*pi*p/ncolss
	
	s <- 1:nrowss
	q <- c( -1*floor((s-1)/2 ), floor(s/2) )
	q <- sort( unique(q) )
	om2 <- 2*pi*q/nrowss
	
	ffs <- expand.grid(om1, om2)
	return(cbind(ffs[,1], ffs[,2]))
}
#' @keywords internal
periodogram = function(spdata)
{
	nrows <- length(unique(mydata[,1])); ncols <- length(unique(mydata[,2]));
	lags <- cov_lags(nrows, ncols)
	spdata <- scale_coords(spdata)
	chat <- est_cov(lags, spdata, nrows, ncols)
	cdata <- cov_complete(chat, lags, ret.mat = F)
	freqs <- get_Fourier_freqs(nrows, ncols)
	nfreqs <- dim(freqs)[1]
	ivec <- c()
	for(i in 1:nfreqs)
	{
		mysum <- 0
		mysum <- sum( cdata[,3]*cos(freqs[i,1]*cdata[,1] + freqs[i,2]*cdata[,2]) )
		I <- mysum#/((2*pi)^2)
		ivec <- c(ivec, I)
	}
	
	return(cbind(freqs, ivec))
}

periodogram2 = function(cdata)
{
	data <- fft(cdata)
	xCx <- data*Conj(data)
	return(xCx/nrow(xCx)/ncol(xCx))#/((2*pi)^2))
}

fftshift <- function(input_matrix, dim = -1) {

    rows <- dim(input_matrix)[1]    
    cols <- dim(input_matrix)[2]    

    swap_up_down <- function(input_matrix) {
        rows_half <- ceiling(rows/2)
        return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
    }

    swap_left_right <- function(input_matrix) {
        cols_half <- ceiling(cols/2)
        return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
    }

    if (dim == -1) {
        input_matrix <- swap_up_down(input_matrix)
        return(swap_left_right(input_matrix))
    }
    else if (dim == 1) {
        return(swap_up_down(input_matrix))
    }
    else if (dim == 2) {
        return(swap_left_right(input_matrix))
    }
    else {
        stop("Invalid dimension parameter")
    }
}

ifftshift <- function(input_matrix, dim = -1) {

    rows <- dim(input_matrix)[1]    
    cols <- dim(input_matrix)[2]    

    swap_up_down <- function(input_matrix) {
        rows_half <- floor(rows/2)
        return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
    }

    swap_left_right <- function(input_matrix) {
        cols_half <- floor(cols/2)
        return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
    }

    if (dim == -1) {
        input_matrix <- swap_left_right(input_matrix)
        return(swap_up_down(input_matrix))
    }
    else if (dim == 1) {
        return(swap_up_down(input_matrix))
    }
    else if (dim == 2) {
        return(swap_left_right(input_matrix))
    }
    else {
        stop("Invalid dimension parameter")
    }
}

periodogram.fn = function(omega, lam, loc.mat, dat.vec, mean.corr = F)
                  {

                   if(mean.corr==T)
                       {
                            dat.vec = dat.vec -mean(dat.vec)
                       }

                    N = length(dat.vec)


                   tmp.site.mat.1 = matrix(rep(loc.mat[,1], N), nrow =N, byrow = T)
                   diff.site.mat.1 = t(tmp.site.mat.1)-tmp.site.mat.1

                   tmp.site.mat.2 = matrix(rep(loc.mat[,2], N), nrow =N, byrow = T)
                   diff.site.mat.2 = t(tmp.site.mat.2)-tmp.site.mat.2

                   tmp.mat = complex(real =cos(diff.site.mat.1*omega[1]+diff.site.mat.2*omega[2]), imaginary=sin(diff.site.mat.1*omega[1]+diff.site.mat.2*omega[2]))

                   periodogram.untap = (dat.vec%*%t(dat.vec))*tmp.mat
                   periodogram = (lam^2)*sum(periodogram.untap)/N#^2


                   periodogram.adj = periodogram-(lam^2)*t(dat.vec%*%dat.vec)/N#^2

                   out = periodogram#.adj
                   return(Re(out))
                   }
