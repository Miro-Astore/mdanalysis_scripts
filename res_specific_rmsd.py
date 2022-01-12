#usage: python mdanalysis/res_specific_rmsd.py TOPFILE TRAJFILE REF_STRUCTURE OUT_NAME.DAT
import numpy as np
import os 
import MDAnalysis as mda
import MDAnalysis.analysis.rms
import sys
#system_list=np.loadtxt("system_list",dtype=str)
cwd=os.getcwd()

#f=open('system_list','r')
#f1=f.readlines()


#os.chdir ('..')
#print (i)
#files=i.split()
#psf=files[0]
#xtc=files[1]
psf=sys.argv[1]
xtc=sys.argv[2]
ref=sys.argv[3]
r = mda.Universe(ref)
out_name=sys.argv[4]
u = mda.Universe(psf,xtc)
R = MDAnalysis.analysis.rms.RMSD(u,r,select=
                                 "name CA and resid 1-70 and not resname UNK",
                                # "name CA and resid 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 856 857 858 859 860 861 862 863 864 865 866 867 868 869 870 871 872 873 874 875 876 877 878 879 880 881 882 883 884 885 886 887 888 889 900 901 902 903 904 905 906 907 908 909 910 911 912 913 914 915 916 917 918 919 920 921 922 923 924 925 926 927 928 929 930 931 932 933 934 935 936 937 938 939 940 941 942 943 944 945 946 947 948 949 950 951 952 953 954 955 956 957 958 959 960 961 962 963 964 965 966 967 968 969 970 971 972 973 974 975 976 977 978 979 980 981 982 983 984 985 986 987 988 989 990 991 992 993 994 995 996 997 998 999 1000 1001 1002 1003 1004 1005 1006 1007 1008 1009 1010 1011 1012 1013 1014 1015 1016 1017 1018 1019 1020 1021 1022 1023 1024 1025 1026 1027 1028 1029 1030 1031 1032 1033 1034 1035 1036 1037 1038 1039 1040 1041 1042 1043 1044 1045 1046 1047 1048 1049 1050 1051 1052 1053 1054 1055 1056 1057 1058 1059 1060 1061 1062 1063 1064 1065 1066 1067 1068 1069 1070 1071 1072 1073 1074 1075 1076 1077 1078 1079 1080 1081 1082 1083 1084 1085 1086 1087 1088 1089 1090 1091 1092 1093 1094 1095 1096 1097 1098 1099 1100 1101 1102 1103 1104 1105 1106 1107 1108 1109 1110 1111 1112 1113 1114 1115 1116 1117 1118 1119 1120 1121 1122 1123 1124 1125 1126 1127 1128 1129 1130 1131 1132 1133 1134 1135 1136 1137 1138 1139 1140 1141 1142 1143 1144 1145 1146 1147 1148 1149 1150 1151 1152 1153 1154 1155 1156 1157 1158 1159 1160 1161 1162 1163 1164 1165 1166 1167 and not resname UNK",
                 groupselections=["name CA and resid 1-44 and not resname UNK ",
                                 "name CA and resid 1 and not resname UNK ",
                                 "name CA and resid 2 and not resname UNK ",
                                 "name CA and resid 3 and not resname UNK ",
                                 "name CA and resid 4 and not resname UNK ",
                                 "name CA and resid 5 and not resname UNK ",
                                 "name CA and resid 6 and not resname UNK ",
                                 "name CA and resid 7 and not resname UNK ",
                                 "name CA and resid 8 and not resname UNK ",
                                 "name CA and resid 9 and not resname UNK ",
                                 "name CA and resid 10 and not resname UNK ",
                                 "name CA and resid 11 and not resname UNK ",
                                 "name CA and resid 12 and not resname UNK ",
                                 "name CA and resid 13 and not resname UNK ",
                                 "name CA and resid 14 and not resname UNK ",
                                 "name CA and resid 15 and not resname UNK ",
                                 "name CA and resid 16 and not resname UNK ",
                                 "name CA and resid 17 and not resname UNK ",
                                 "name CA and resid 18 and not resname UNK ",
                                 "name CA and resid 19 and not resname UNK ",
                                 "name CA and resid 20 and not resname UNK ",
                                 "name CA and resid 21 and not resname UNK ",
                                 "name CA and resid 22 and not resname UNK ",
                                 "name CA and resid 23 and not resname UNK ",
                                 "name CA and resid 24 and not resname UNK ",
                                 "name CA and resid 25 and not resname UNK ",
                                 "name CA and resid 26 and not resname UNK ",
                                 "name CA and resid 27 and not resname UNK ",
                                 "name CA and resid 28 and not resname UNK ",
                                 "name CA and resid 29 and not resname UNK ",
                                 "name CA and resid 30 and not resname UNK ",
                                 "name CA and resid 31 and not resname UNK ",
                                 "name CA and resid 32 and not resname UNK ",
                                 "name CA and resid 33 and not resname UNK ",
                                 "name CA and resid 34 and not resname UNK ",
                                 "name CA and resid 35 and not resname UNK ",
                                 "name CA and resid 36 and not resname UNK ",
                                 "name CA and resid 37 and not resname UNK ",
                                 "name CA and resid 38 and not resname UNK ",
                                 "name CA and resid 39 and not resname UNK ",
                                 "name CA and resid 40 and not resname UNK ",
                                 "name CA and resid 41 and not resname UNK ",
                                 "name CA and resid 42 and not resname UNK ",
                                 "name CA and resid 43 and not resname UNK ",
                                 "name CA and resid 44 and not resname UNK ",
                                 "name CA and resid 45 and not resname UNK ",
                                 "name CA and resid 46 and not resname UNK ",
                                 "name CA and resid 47 and not resname UNK ",
                                 "name CA and resid 48 and not resname UNK ",
                                 "name CA and resid 49 and not resname UNK ",
                                 "name CA and resid 50 and not resname UNK ",
                                 "name CA and resid 51 and not resname UNK ",
                                 "name CA and resid 52 and not resname UNK ",
                                 "name CA and resid 53 and not resname UNK ",
                                 "name CA and resid 54 and not resname UNK ",
                                 "name CA and resid 55 and not resname UNK ",
                                 "name CA and resid 56 and not resname UNK ",
                                 "name CA and resid 57 and not resname UNK ",
                                 "name CA and resid 58 and not resname UNK ",
                                 "name CA and resid 59 and not resname UNK ",
                                 "name CA and resid 60 and not resname UNK ",
                                 "name CA and resid 61 and not resname UNK ",
                                 "name CA and resid 62 and not resname UNK ",
                                 "name CA and resid 63 and not resname UNK ",
                                 "name CA and resid 64 and not resname UNK ",
                                 "name CA and resid 65 and not resname UNK ",
                                 "name CA and resid 66 and not resname UNK ",
                                 "name CA and resid 67 and not resname UNK ",
                                 "name CA and resid 68 and not resname UNK ",
                                 "name CA and resid 69 and not resname UNK ",
                                 "name CA and resid 70 and not resname UNK "], filename=out_name)
R.run()
#R.save()                      
np.savetxt(out_name, R.rmsd)
