#!/bin/bash

GREP_OPTIONS=''

cookiejar=$(mktemp cookies.XXXXXXXXXX)
netrc=$(mktemp netrc.XXXXXXXXXX)
chmod 0600 "$cookiejar" "$netrc"
function finish {
  rm -rf "$cookiejar" "$netrc"
}

trap finish EXIT
WGETRC="$wgetrc"

prompt_credentials() {
    echo "Enter your Earthdata Login or other provider supplied credentials"
    read -p "Username (yinonbaron): " username
    username=${username:-yinonbaron}
    read -s -p "Password: " password
    echo "machine urs.earthdata.nasa.gov login $username password $password" >> $netrc
    echo
}

exit_with_error() {
    echo
    echo "Unable to Retrieve Data"
    echo
    echo $1
    echo
    echo "https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1982001_001_2018224204211/VCF5KYR_1982001_001_2018224204211.tif"
    echo
    exit 1
}

prompt_credentials
  detect_app_approval() {
    approved=`curl -s -b "$cookiejar" -c "$cookiejar" -L --max-redirs 5 --netrc-file "$netrc" https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1982001_001_2018224204211/VCF5KYR_1982001_001_2018224204211.tif -w '\n%{http_code}' | tail  -1`
    if [ "$approved" -ne "200" ] && [ "$approved" -ne "301" ] && [ "$approved" -ne "302" ]; then
        # User didn't approve the app. Direct users to approve the app in URS
        exit_with_error "Please ensure that you have authorized the remote application by visiting the link below "
    fi
}

setup_auth_curl() {
    # Firstly, check if it require URS authentication
    status=$(curl -s -z "$(date)" -w '\n%{http_code}' https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1982001_001_2018224204211/VCF5KYR_1982001_001_2018224204211.tif | tail -1)
    if [[ "$status" -ne "200" && "$status" -ne "304" ]]; then
        # URS authentication is required. Now further check if the application/remote service is approved.
        detect_app_approval
    fi
}

setup_auth_wget() {
    # The safest way to auth via curl is netrc. Note: there's no checking or feedback
    # if login is unsuccessful
    touch ~/.netrc
    chmod 0600 ~/.netrc
    credentials=$(grep 'machine urs.earthdata.nasa.gov' ~/.netrc)
    if [ -z "$credentials" ]; then
        cat "$netrc" >> ~/.netrc
    fi
}

fetch_urls() {
  if command -v curl >/dev/null 2>&1; then
      setup_auth_curl
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        curl -f -b "$cookiejar" -c "$cookiejar" -L --netrc-file "$netrc" -g -o $stripped_query_params -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  elif command -v wget >/dev/null 2>&1; then
      # We can't use wget to poke provider server to get info whether or not URS was integrated without download at least one of the files.
      echo
      echo "WARNING: Can't find curl, use wget instead."
      echo "WARNING: Script may not correctly identify Earthdata Login integrations."
      echo
      setup_auth_wget
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        wget --load-cookies "$cookiejar" --save-cookies "$cookiejar" --output-document $stripped_query_params --keep-session-cookies -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  else
      exit_with_error "Error: Could not find a command-line downloader.  Please install curl or wget"
  fi
}

fetch_urls <<'EDSCEOF'
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1982001_001_2018224204211/VCF5KYR_1982001_001_2018224204211.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1983001_001_2018224204436/VCF5KYR_1983001_001_2018224204436.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1985001_001_2018224204542/VCF5KYR_1985001_001_2018224204542.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1984001_001_2018224204542/VCF5KYR_1984001_001_2018224204542.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1986001_001_2018224204727/VCF5KYR_1986001_001_2018224204727.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1987001_001_2018224204756/VCF5KYR_1987001_001_2018224204756.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1988001_001_2018224204838/VCF5KYR_1988001_001_2018224204838.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1989001_001_2018224204904/VCF5KYR_1989001_001_2018224204904.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1990001_001_2018224204938/VCF5KYR_1990001_001_2018224204938.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1991001_001_2018224205008/VCF5KYR_1991001_001_2018224205008.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1992001_001_2018224205034/VCF5KYR_1992001_001_2018224205034.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1993001_001_2018224205123/VCF5KYR_1993001_001_2018224205123.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1995001_001_2018224205251/VCF5KYR_1995001_001_2018224205251.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1996001_001_2018224205332/VCF5KYR_1996001_001_2018224205332.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1997001_001_2018224205415/VCF5KYR_1997001_001_2018224205415.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1998001_001_2018224205440/VCF5KYR_1998001_001_2018224205440.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_1999001_001_2018224205513/VCF5KYR_1999001_001_2018224205513.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2001001_001_2018224205557/VCF5KYR_2001001_001_2018224205557.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2002001_001_2018224205629/VCF5KYR_2002001_001_2018224205629.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2003001_001_2018224205710/VCF5KYR_2003001_001_2018224205710.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2004001_001_2018224205736/VCF5KYR_2004001_001_2018224205736.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2005001_001_2018224205812/VCF5KYR_2005001_001_2018224205812.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2006001_001_2018224205841/VCF5KYR_2006001_001_2018224205841.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2007001_001_2018224205909/VCF5KYR_2007001_001_2018224205909.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2008001_001_2018224205932/VCF5KYR_2008001_001_2018224205932.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2009001_001_2018224210002/VCF5KYR_2009001_001_2018224210002.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2010001_001_2018224210033/VCF5KYR_2010001_001_2018224210033.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2011001_001_2018224210105/VCF5KYR_2011001_001_2018224210105.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2012001_001_2018224210128/VCF5KYR_2012001_001_2018224210128.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2013001_001_2018224210156/VCF5KYR_2013001_001_2018224210156.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2014001_001_2018224210222/VCF5KYR_2014001_001_2018224210222.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2015001_001_2018224210248/VCF5KYR_2015001_001_2018224210248.tif
https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/VCF5KYR.001/VCF5KYR_2016001_001_2018224210310/VCF5KYR_2016001_001_2018224210310.tif
EDSCEOF