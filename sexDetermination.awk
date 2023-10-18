#!/bin/bash
###########################################################################################
###########################################################################################
#
# Author : Patrick Cichowicz
# 
# Masters Thesis project script for aDNA and UDG (USER) enzyme reduction treatment analysis
# University of Copenhagen, GLOBE Institute, SUND
# 
###########################################################################################
###########################################################################################

BEGIN {
    xReads = 0
    yReads = 0
    mReads = 0
    autReads = 0

    xSites = 0
    ySites = 0
    mSites = 0
    autSites = 0

    xCov = 0
    yCov = 0
    mCov = 0
    autCov = 0
}
{
    chr = $1
    pos = $2
    cov = $3

    if(chr == "chrX") {
        xReads += cov
        xSites += 1 
        if (cov > 0) xCov +=1
    }
    else if(chr == "chrY") {
        yReads += cov
        ySites += 1
        if (cov > 0) yCov +=1
    }
    else if(chr == "chrM") {
        mReads += cov
        mSites += 1
        if (cov > 0) mCov +=1
    }
    else {
        autReads += cov
        autSites += 1
        if (cov > 0) autCov +=1
    }
}
END {
    OFS="\t"
    print("Chr" ,"Sites", "Depth")
    print("xDepth", xCov, xSites > 0 ? xReads / xSites : 0)
    print("ydepth", yCov, ySites > 0 ? yReads / ySites : 0)
    print("mDepth", mCov, mSites > 0 ? mReads / mSites : 0)
    print("autDepth", autCov, autSites > 0 ? autReads / autSites : 0)
}
