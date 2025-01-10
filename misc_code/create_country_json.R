
## Prepare vocabulary for PhyloTwin frontend
# `client/src/Vocabularies/country.json`

## Country names should match polygons from Natural Earth v.5.1.2 [May 13, 2022]
# https://www.naturalearthdata.com/downloads/

supported_codes <- c("AD", "AE", "AF", "AG", "AI", "AL", "AM", "AO", "AQ", "AR", 
	"AS", "AT", "AU", "AW", "AX", "AZ", "BA", "BB", "BD", "BE", "BF", 
	"BG", "BH", "BI", "BJ", "BL", "BM", "BN", "BO", "BR", "BS", "BT", 
	"BW", "BY", "BZ", "CA", "CD", "CF", "CG", "CH", "CI", "CK", "CL", 
	"CM", "CN", "CO", "CR", "CU", "CV", "CW", "CY", "CZ", "DE", "DJ", 
	"DK", "DM", "DO", "DZ", "EC", "EE", "EG", "EH", "ER", "ES", "ET", 
	"FI", "FJ", "FK", "FM", "FO", "FR", "GA", "GB", "GD", "GE", "GG", 
	"GH", "GL", "GM", "GN", "GQ", "GR", "GS", "GT", "GU", "GW", "GY", 
	"HK", "HM", "HN", "HR", "HT", "HU", "ID", "IE", "IL", "IM", "IN", 
	"IO", "IQ", "IR", "IS", "IT", "JE", "JM", "JO", "JP", "KE", "KG", 
	"KH", "KI", "KM", "KN", "KP", "KR", "KW", "KY", "KZ", "LA", "LB", 
	"LC", "LI", "LK", "LR", "LS", "LT", "LU", "LV", "LY", "MA", "MC", 
	"MD", "ME", "MF", "MG", "MH", "MK", "ML", "MM", "MN", "MO", "MP", 
	"MR", "MS", "MT", "MU", "MV", "MW", "MX", "MY", "MZ", "NA", "NC", 
	"NE", "NF", "NG", "NI", "NL", "NO", "NP", "NR", "NU", "NZ", "OM", 
	"PA", "PE", "PF", "PG", "PH", "PK", "PL", "PM", "PN", "PR", "PS", 
	"PT", "PW", "PY", "QA", "RO", "RS", "RU", "RW", "SA", "SB", "SC", 
	"SD", "SE", "SG", "SH", "SI", "SK", "SL", "SM", "SN", "SO", "SR", 
	"SS", "ST", "SV", "SX", "SY", "SZ", "TC", "TD", "TF", "TG", "TH", 
	"TJ", "TL", "TM", "TN", "TO", "TR", "TT", "TV", "TW", "TZ", "UA", 
	"UG", "US", "UY", "UZ", "VA", "VC", "VE", "VG", "VI", "VN", "VU", 
	"WF", "WS", "XK", "YE", "ZA", "ZM", "ZW")


library(data.table)
library(sf)
library(jsonlite)

## Load the Natural Earth data
# cd /tmp
# wget https://naciscdn.org/naturalearth/packages/natural_earth_vector.gpkg.zip
# unzip natural_earth_vector.gpkg.zip
NE <- st_read("/tmp/packages/natural_earth_vector.gpkg", layer = "ne_50m_admin_0_countries")
NES <- NE[ ! NE$ISO_A2_EH %in% "-99", c("NAME", "ISO_A2_EH") ]
setDT(NES)
NES[, geom := NULL ]
NES <- unique(NES)
setnames(NES, new = c("CountryUNICODE", "A2"))

if(any(! supported_codes %in% NES$A2)){
  supported_codes[ ! supported_codes %in% NES$A2 ]
}

NES[, Country := iconv(CountryUNICODE, "UTF-8", "ASCII//TRANSLIT") ]
NES[  Country %in% "eSwatini" , Country := "Eswatini" ]

setorder(NES, Country)


json_data <- toJSON(NES[ , .(Country, A2)], pretty = TRUE)
write(json_data, file = "country.json")



############################################### misc

library(data.table)

# https://github.com/nocworx/datasets-country-codes/
# https://www.statoids.com/wab.html

tab <- fread("https://raw.githubusercontent.com/nocworx/datasets-country-codes/refs/heads/main/data/country-codes.csv")

if(any(! supported_codes %in% tab$"ISO3166-1-Alpha-2")){
  supported_codes[ ! supported_codes %in% tab$"ISO3166-1-Alpha-2" ]     # "NA" "XK"
}

tmp <- tab[ , .(`official_name_en`, `ISO3166-1-Alpha-2`) ]
setnames(x = tmp, new = c("Country", "A2"))

res <- tmp[ A2 %in% supported_codes ]

add <- rowwiseDT(
  Country=, A2=,
  "Namibia" ,"NA",
  "Kosovo", "XK")

res <- rbind(res, add)


