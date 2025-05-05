import pandas as pd
import numpy as np

ds = pd.read_csv('.csv file path')

a = 6378137
b = 6356752.3142
esq = 6.69437999014 * 0.001
e1sq = 6.73949674228 * 0.001

def ecef2geodetic(ecef, radians=False):
  """
  Convert ECEF coordinates to geodetic using ferrari's method
  """
  # Save shape and export column
  ecef = np.atleast_1d(ecef)
  input_shape = ecef.shape
  ecef = np.atleast_2d(ecef)
  x, y, z = ecef[:, 0], ecef[:, 1], ecef[:, 2]

  ratio = 1.0 if radians else (180.0 / np.pi)

  # Conver from ECEF to geodetic using Ferrari's methods
  # https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#Ferrari.27s_solution
  r = np.sqrt(x * x + y * y)
  Esq = a * a - b * b
  F = 54 * b * b * z * z
  G = r * r + (1 - esq) * z * z - esq * Esq
  C = (esq * esq * F * r * r) / (pow(G, 3))
  S = np.cbrt(1 + C + np.sqrt(C * C + 2 * C))
  P = F / (3 * pow((S + 1 / S + 1), 2) * G * G)
  Q = np.sqrt(1 + 2 * esq * esq * P)
  r_0 =  -(P * esq * r) / (1 + Q) + np.sqrt(0.5 * a * a*(1 + 1.0 / Q) - \
        P * (1 - esq) * z * z / (Q * (1 + Q)) - 0.5 * P * r * r)
  U = np.sqrt(pow((r - esq * r_0), 2) + z * z)
  V = np.sqrt(pow((r - esq * r_0), 2) + (1 - esq) * z * z)
  Z_0 = b * b * z / (a * V)
  h = U * (1 - b * b / (a * V))
  lat = ratio*np.arctan((z + e1sq * Z_0) / r)
  lon = ratio*np.arctan2(y, x)

  # stack the new columns and return to the original shape
  geodetic = np.column_stack((lon, lat, h))
  return geodetic.reshape(input_shape)


def generate_kml(group_data, prn):
    kml_content = [
        '<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n<Document>\n<name>Paths_{prn}</name>\n'
    ]
    
    kml_content.extend([
        "<Style id=\"0\">\n<LineStyle>\n<color>ff00fd00</color>\n<width>1.</width>\n</LineStyle>\n</Style>\n",
        "<Style id=\"1\">\n<LineStyle>\n<color>ff0000ff</color>\n<width>1.</width>\n</LineStyle>\n</Style>\n",
        "<Style id=\"2\">\n<LineStyle>\n<color>ff0083ff</color>\n<width>1.</width>\n</LineStyle>\n</Style>\n"
    ])
    
    for index, row in group_data.iterrows():
        receivergt = row[['GT_Longitude', 'GT_Latitude', 'GT_Height']].values
        sd = row[['SaX', 'SaY', 'SaZ']].values
        sv = ecef2geodetic(sd)
        #接收机位置
        kml_content.append(f"""<Placemark>
          <name>Receiver {index}</name>
          <Point>
            <altitudeMode>absolute</altitudeMode>
            <coordinates>{receivergt[0]}, {receivergt[1]}, {receivergt[2]}</coordinates>
          </Point>
        </Placemark>""")

        # 添加路径
        outs = f"{receivergt[0]},{receivergt[1]},{receivergt[2]} {sv[0]},{sv[1]},{sv[2]}"
        kml_content.append(f'<Placemark>\n<open>1</open>\n<styleUrl>#2</styleUrl>\n<LineString>\n<tessellate>1</tessellate>\n<altitudeMode>absolute</altitudeMode>\n<coordinates>\n{outs}\n</coordinates>\n</LineString>\n</Placemark>\n')
    
    kml_content.extend([
        '</Document>\n</kml>'
    ])
    
    return ''.join(kml_content)
groups = ds.groupby('PRN')

# 指定保存KML文件的文件夹路径
output_folder = "C:/Users/ZHANG/OneDrive - The Hong Kong Polytechnic University/桌面/NewData/New_Path/1side_diff_d2a/"

# 遍历每个PRN组并生成KML文件
for prn, group_data in groups:
    kml_content = generate_kml(group_data, prn)
    file_name = f'{output_folder}receiver_satellite_paths_PRN_{prn}.kml'
    with open(file_name, 'w') as f:
        f.write(kml_content)
    print(f"KML文件已生成并保存为{file_name}")
