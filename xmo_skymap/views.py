import io
import os

import matplotlib
import numpy as np
import pytz

from astropy.coordinates import EarthLocation, AltAz, get_sun, Angle, get_body, SkyCoord
from astropy.time import Time
from astropy.visualization import astropy_mpl_style, quantity_support
from django.http import HttpResponse, JsonResponse
from datetime import datetime, timedelta
from django.core.cache import cache
from django.template import loader
import astropy.units as u

from xmo_skymap.plan_gen_tool import load_item_template, write_plan_file
from xmo_skymap.rotate_tools import get_l_r_t_b_axis, rotate, get_top, get_left, get_right, get_bottom, \
    get_rotate_fix_axis, get_top_fix_axis, get_right_fix_axis, get_bottom_fix_axis, get_left_fix_axis

matplotlib.use('TkAgg')
from matplotlib import pyplot as plt


def index(request):
    # return HttpResponse("skymap")
    now_date = datetime.today().strftime('%Y%m%d')
    req_date = request.GET.get('req_date', now_date)
    print(' -> req [%s]' % req_date)

    latest_tele_system_list = []

    print(cache.get(req_date))
    sky_process_map = cache.get(req_date)
    template = loader.get_template('skymap/aladin.html')
    context = {
        'latest_tele_system_list': latest_tele_system_list,
        # 'day_in_cache': num,
        'process_map': sky_process_map,
        'req_date': req_date,

    }
    return HttpResponse(template.render(context, request))


def draw_sun(request):
    now_date = datetime.today().strftime('%Y%m%d')
    url = request.GET.get('url_for_test', now_date)
    req_overlap = float(request.GET.get('req_overlap_deg', 0.2))
    req_row = int(request.GET.get('req_row_int', 2))
    req_col = int(request.GET.get('req_col_int', 2))
    req_x_img = float(request.GET.get('req_x_img_deg', 2.1))
    req_y_img = float(request.GET.get('req_y_img_deg', 4.2))
    req_ra_cen_deg = float(request.GET.get('req_ra_deg', 30.1))
    req_dec_cen_deg = float(request.GET.get('req_dec_deg', 40.3))
    ex_message = ''

    # 设置观测地点（以纽约为例）
    location = EarthLocation(lat='39.7128', lon='116.0060')  # 纬度和经度

    # 获取当前时间
    current_time = datetime.now()
    # 获取次日凌晨0点0分的时间
    current_time = current_time.replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(days=1)
    next_day_midnight = Time(current_time)

    # 获取当前时间
    current_time = Time.now()
    print(f"time {next_day_midnight.iso}")
    print(f"当前时间：{next_day_midnight.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    current_time_with_tz = current_time.to_datetime(pytz.timezone('Asia/Shanghai'))
    # 将观测地点和时间结合起来
    altaz = AltAz(location=location, obstime=current_time)
    print(f"当前时间：{current_time_with_tz.strftime('%Y-%m-%d %H:%M:%S %Z')}")

    # 生成前后各12小时的时间数组
    time_array = []
    sun_curve = []
    moon_curve = []
    for i in range(-12, 13):
        time_offset = i * u.hour
        new_time = next_day_midnight + time_offset
        time_array.append(new_time)
        # print(new_time)
        line_sun = get_body('sun', new_time, location)
        line_sun_radec = line_sun.transform_to(AltAz(location=location, obstime=new_time)).transform_to('gcrs')
        line_moon = get_body('moon', new_time, location)
        line_moon_radec = line_moon.transform_to(AltAz(location=location, obstime=new_time)).transform_to('gcrs')
        sun_curve.append([str(line_sun_radec.ra.deg), str(line_sun_radec.dec.deg)])
        moon_curve.append([str(line_moon_radec.ra.deg), str(line_moon_radec.dec.deg)])

    # 计算太阳的赤经和赤纬
    sun = get_body('sun', next_day_midnight, location)
    print(sun)
    sun_radec = sun.transform_to(AltAz(location=location, obstime=next_day_midnight)).transform_to('gcrs')

    # 计算月亮的赤经和赤纬
    moon = get_body('moon', next_day_midnight, location)
    moon_radec = moon.transform_to(AltAz(location=location, obstime=next_day_midnight)).transform_to('gcrs')
    print(moon)
    sun_ra_hms = Angle(sun_radec.ra.deg, unit=u.deg)
    moon_ra_hms = Angle(moon_radec.ra.deg, unit=u.deg)

    # ra_hms_str = ra_hms.to_string(unit=u.hour, sep=':', precision=2)
    # 打印太阳和月亮的赤经和赤纬
    print(f"太阳的赤经（Right Ascension）：{sun_radec.ra.deg}  {sun_ra_hms.hms}  {sun_ra_hms.dms}")
    print(f"太阳的赤纬（Declination）：{sun_radec.dec.deg}")
    print()
    print(f"月亮的赤经（Right Ascension）：{moon_radec.ra.deg}   {moon_ra_hms.hms}   {moon_ra_hms.dms}")
    print(f"月亮的赤纬（Declination）：{moon_radec.dec.deg}")
    print(f"time {current_time.iso}")
    print(f"当前时间：{next_day_midnight.strftime('%Y-%m-%d %H:%M:%S %Z')}")

    # 起始中心坐标（示例坐标，你可以根据实际需要修改）
    center_ra = req_ra_cen_deg * u.deg
    center_dec = req_dec_cen_deg * u.deg

    square_list = []
    center_list = []

    # 5x4 图像，间距重叠1 中心上下间距7 左右间距9
    num_row = req_row
    num_colum = req_col
    img_wid = req_x_img
    img_hei = req_y_img
    img_overlap = req_overlap
    row_head_coord_list = []
    cen_ra_row = center_ra
    cen_dec_row = center_dec
    coord_debug = SkyCoord(ra=center_ra, dec=center_dec, unit='deg')
    az_alt = coord_debug.transform_to(AltAz(obstime=current_time, location=location))
    print('=========   az  [%s]    alt [%s]   ========' % (az_alt.az, az_alt.alt))
    cart_debug = coord_debug.cartesian
    print('========= xyz  (%s,%s,%s)   ========' % (cart_debug.x.value, cart_debug.y.value, cart_debug.z.value))

    # =======================================================================================
    rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(cen_ra_row, cen_dec_row)
    rtt_l = np.array([0, 0, 1])
    rtt_r = np.array([0, 0, -1])
    # rtt_t = np.array([1, 0, 0])
    # rtt_b = np.array([-1, 0, 0])
    if num_row % 2 == 0:
        half_row = num_row // 2
        for i in range(half_row):
            cord_row_head_center_item = get_top_fix_axis(cen_ra_row, cen_dec_row, ((img_hei - img_overlap) * i)
                                                         + ((img_hei / 2) - (img_overlap / 2)), rtt_t)
            row_head_coord_list.append(cord_row_head_center_item)
        # reverse top list
        row_head_coord_list.reverse()
        for i in range(half_row):
            cord_row_head_center_item = get_bottom_fix_axis(cen_ra_row, cen_dec_row, ((img_hei - img_overlap) * i)
                                                            + ((img_hei / 2) - (img_overlap / 2)), rtt_b)
            row_head_coord_list.append(cord_row_head_center_item)
    else:
        half_row = num_row // 2
        for i in range(half_row):
            cord_row_head_center_item = get_top_fix_axis(cen_ra_row, cen_dec_row, (img_hei - img_overlap) * (i + 1),
                                                         rtt_t)
            row_head_coord_list.append(cord_row_head_center_item)
        # reverse top list
        row_head_coord_list.reverse()
        row_head_coord_list.append(SkyCoord(ra=cen_ra_row, dec=cen_dec_row, unit='deg').transform_to('gcrs'))
        for i in range(half_row):
            cord_row_head_center_item = get_bottom_fix_axis(cen_ra_row, cen_dec_row, (img_hei - img_overlap) * (i + 1),
                                                            rtt_b)
            row_head_coord_list.append(cord_row_head_center_item)

    # 一种起始点为左下角的row算法
    # for i in range(num_row):
    #     cord_row_head_center_item = get_top_fix_axis(cen_ra_row, cen_dec_row, (img_hei-img_overlap)*i, rtt_t)
    #     row_head_coord_list.append(cord_row_head_center_item)
    # 计算所有实际中心点

    for i in range(num_row):
        row_center_ra = row_head_coord_list[i].ra.value
        row_center_dec = row_head_coord_list[i].dec.value
        cord_t_c = row_head_coord_list[i]

        if num_colum % 2 == 0:
            half_colum = num_colum // 2
            temp_center_coordinates = []
            pre_center_ra_colum = row_center_ra
            pre_center_dec_colum = row_center_dec
            for j in range(half_colum):
                # rtt_l = np.array([0, 0, 1])
                # rtt_r = np.array([0, 0, -1])
                rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(pre_center_ra_colum, pre_center_dec_colum)
                first_angle = 1
                if j == 0:
                    first_angle = 2
                cord_t_c_dec_by_z = get_left_fix_axis(pre_center_ra_colum, pre_center_dec_colum,
                                                      ((img_wid - img_overlap) / first_angle), rtt_l)
                pre_center_ra_colum = cord_t_c_dec_by_z.ra.value
                # pre_center_dec_colum = cord_t_c_dec_by_z.dec.value
                temp_center_coordinates.append([pre_center_ra_colum, pre_center_dec_colum])
            temp_center_coordinates.reverse()
            pre_center_ra_colum = row_center_ra
            pre_center_dec_colum = row_center_dec
            for j in range(half_colum):
                # rtt_l = np.array([0, 0, 1])
                # rtt_r = np.array([0, 0, -1])
                rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(pre_center_ra_colum, pre_center_dec_colum)
                first_angle = 1
                if j == 0:
                    first_angle = 2
                cord_t_c_dec_by_z = get_right_fix_axis(pre_center_ra_colum, pre_center_dec_colum,
                                                       ((img_wid - img_overlap) / first_angle), rtt_r)
                pre_center_ra_colum = cord_t_c_dec_by_z.ra.value
                # pre_center_dec_colum = cord_t_c_dec_by_z.dec.value
                temp_center_coordinates.append([pre_center_ra_colum, pre_center_dec_colum])
            center_list.extend(temp_center_coordinates)
        else:
            half_colum = num_colum // 2
            temp_center_coordinates = []
            pre_center_ra_colum = row_center_ra
            pre_center_dec_colum = row_center_dec
            for j in range(half_colum):
                # rtt_l = np.array([0, 0, 1])
                # rtt_r = np.array([0, 0, -1])
                rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(pre_center_ra_colum, pre_center_dec_colum)
                cord_t_c_dec_by_z = get_left_fix_axis(pre_center_ra_colum, pre_center_dec_colum,
                                                      (img_wid - img_overlap), rtt_l)
                pre_center_ra_colum = cord_t_c_dec_by_z.ra.value
                # pre_center_dec_colum = cord_t_c_dec_by_z.dec.value
                temp_center_coordinates.append([pre_center_ra_colum, pre_center_dec_colum])
            temp_center_coordinates.reverse()
            temp_center_coordinates.append([row_center_ra, row_center_dec])
            pre_center_ra_colum = row_center_ra
            pre_center_dec_colum = row_center_dec
            for j in range(half_colum):
                # rtt_l = np.array([0, 0, 1])
                # rtt_r = np.array([0, 0, -1])
                rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(pre_center_ra_colum, pre_center_dec_colum)
                cord_t_c_dec_by_z = get_right_fix_axis(pre_center_ra_colum, pre_center_dec_colum,
                                                       (img_wid - img_overlap), rtt_r)
                pre_center_ra_colum = cord_t_c_dec_by_z.ra.value
                # pre_center_dec_colum = cord_t_c_dec_by_z.dec.value
                temp_center_coordinates.append([pre_center_ra_colum, pre_center_dec_colum])
            center_list.extend(temp_center_coordinates)

    # 计算四角点
    for i in range(len(center_list)):
        cord_real_center = center_list[i]
        rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(cord_real_center[0], cord_real_center[1])
        # rtt_l = np.array([0, 0, 1])
        # rtt_r = np.array([0, 0, -1])
        # print('[%s]  [%s]              [%s]   [%s]' % (i, j, center_ra_colum, center_dec_colum))
        coordinates = []
        cord_t_c_tm = get_top_fix_axis(cord_real_center[0], cord_real_center[1], (img_hei / 2), rtt_t)
        cord_t_c_bm = get_bottom_fix_axis(cord_real_center[0], cord_real_center[1], (img_hei / 2), rtt_b)
        cord_t_c_tl = get_left_fix_axis(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, (img_wid / 2), rtt_l)
        cord_t_c_tr = get_right_fix_axis(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, (img_wid / 2), rtt_r)
        cord_t_c_bl = get_left_fix_axis(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, (img_wid / 2), rtt_l)
        cord_t_c_br = get_right_fix_axis(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, (img_wid / 2), rtt_r)
        coordinates.append([cord_t_c_tl.ra.value, cord_t_c_tl.dec.value])
        coordinates.append([cord_t_c_tr.ra.value, cord_t_c_tr.dec.value])
        coordinates.append([cord_t_c_br.ra.value, cord_t_c_br.dec.value])
        coordinates.append([cord_t_c_bl.ra.value, cord_t_c_bl.dec.value])

        square_list.append(coordinates)

    # for i in range(num_row):
    #     center_ra_colum = row_head_coord_list[i].ra.value
    #     center_dec_colum = row_head_coord_list[i].dec.value
    #     rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(center_ra_colum, center_dec_colum)
    #     # rtt_l = np.array([0, 0, 1])
    #     # rtt_r = np.array([0, 0, -1])
    #     cord_t_c = row_head_coord_list[i]
    #     for j in range(num_colum):
    #         # cord_t_c = get_right_fix_axis(center_ra_colum, center_dec_colum, (img_wid-img_overlap) * j, rtt_r)
    #         # cord_t_c = get_right_fix_axis(center_ra_colum, center_dec_colum, (img_wid-img_overlap), rtt_r)
    #
    #         rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(cord_t_c.ra.value, cord_t_c.dec.value)
    #         # rtt_l = np.array([0, 0, 1])
    #         # rtt_r = np.array([0, 0, -1])
    #         # print('[%s]  [%s]              [%s]   [%s]' % (i, j, center_ra_colum, center_dec_colum))
    #         coordinates = []
    #         cord_t_c_tm = get_top_fix_axis(cord_t_c.ra.value, cord_t_c.dec.value, (img_hei / 2), rtt_t)
    #         cord_t_c_bm = get_bottom_fix_axis(cord_t_c.ra.value, cord_t_c.dec.value, (img_hei / 2), rtt_b)
    #         cord_t_c_tl = get_left_fix_axis(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, (img_wid / 2), rtt_l)
    #         cord_t_c_tr = get_right_fix_axis(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, (img_wid / 2), rtt_r)
    #         cord_t_c_bl = get_left_fix_axis(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, (img_wid / 2), rtt_l)
    #         cord_t_c_br = get_right_fix_axis(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, (img_wid / 2), rtt_r)
    #         coordinates.append([cord_t_c_tl.ra.value, cord_t_c_tl.dec.value])
    #         coordinates.append([cord_t_c_tr.ra.value, cord_t_c_tr.dec.value])
    #         coordinates.append([cord_t_c_br.ra.value, cord_t_c_br.dec.value])
    #         coordinates.append([cord_t_c_bl.ra.value, cord_t_c_bl.dec.value])
    #
    #         square_list.append(coordinates)
    #         center_list.append([cord_t_c.ra.value, cord_t_c.dec.value])
    #
    #         cord_t_c = get_right_fix_axis(center_ra_colum, center_dec_colum, (img_wid - img_overlap), rtt_r)
    #         rtt_l = np.array([0, 0, 1])
    #         rtt_r = np.array([0, 0, -1])
    #         cord_t_c_dec_by_z = get_right_fix_axis(center_ra_colum, center_dec_colum, (img_wid - img_overlap), rtt_r)
    #         center_ra_colum = cord_t_c.ra.value
    #         center_dec_colum = cord_t_c_dec_by_z.dec.value
    # =======================================================================================
    # for i in range(num_row):
    #     cord_row_head_center_item = get_top(cen_ra_row, cen_dec_row, (img_hei-img_overlap)*i, current_time, location)
    #     row_head_coord_list.append(cord_row_head_center_item)
    # for i in range(num_row):
    #     center_ra_colum = row_head_coord_list[i].ra.value
    #     center_dec_colum = row_head_coord_list[i].dec.value
    #     for j in range(num_colum):
    #         cord_t_c = get_right(center_ra_colum, center_dec_colum, (img_wid-img_overlap) * j, current_time, location)
    #         # print('[%s]  [%s]              [%s]   [%s]' % (i, j, center_ra_colum, center_dec_colum))
    #         coordinates = []
    #         cord_t_c_tm = get_top(cord_t_c.ra.value, cord_t_c.dec.value, (img_hei/2), current_time, location)
    #         cord_t_c_bm = get_bottom(cord_t_c.ra.value, cord_t_c.dec.value, (img_hei/2), current_time, location)
    #         cord_t_c_tl = get_left(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, (img_wid/2), current_time, location)
    #         cord_t_c_tr = get_right(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, (img_wid/2), current_time, location)
    #         cord_t_c_bl = get_left(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, (img_wid/2), current_time, location)
    #         cord_t_c_br = get_right(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, (img_wid/2), current_time, location)
    #         coordinates.append([cord_t_c_tl.ra.value, cord_t_c_tl.dec.value])
    #         coordinates.append([cord_t_c_tr.ra.value, cord_t_c_tr.dec.value])
    #         coordinates.append([cord_t_c_br.ra.value, cord_t_c_br.dec.value])
    #         coordinates.append([cord_t_c_bl.ra.value, cord_t_c_bl.dec.value])
    #
    #         square_list.append(coordinates)

    # =======================================================================================
    # for i in range(2):
    #     cord_row_head_center = get_top(center_ra, center_dec, 7*i, current_time, location)
    #     center_ra = cord_t_c.ra.value
    #     center_dec = cord_t_c.dec.value
    #     for j in range(3):
    #         cord_t_c = get_right(center_ra, center_dec, 9*j, current_time, location)
    #         center_ra = cord_t_c.ra.value
    #         center_dec = cord_t_c.dec.value
    #         print('[%s]  [%s]              [%s]   [%s]' % (i, j, center_ra, center_dec))
    #         coordinates = []
    #         cord_t_c_tm = get_top(cord_t_c.ra.value, cord_t_c.dec.value, 2, current_time, location)
    #         cord_t_c_bm = get_bottom(cord_t_c.ra.value, cord_t_c.dec.value, 2, current_time, location)
    #         cord_t_c_tl = get_left(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, 2.5, current_time, location)
    #         cord_t_c_tr = get_right(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, 2.5, current_time, location)
    #         cord_t_c_bl = get_left(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, 2.5, current_time, location)
    #         cord_t_c_br = get_right(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, 2.5, current_time, location)
    #         coordinates.append([cord_t_c_tl.ra.value, cord_t_c_tl.dec.value])
    #         coordinates.append([cord_t_c_tr.ra.value, cord_t_c_tr.dec.value])
    #         coordinates.append([cord_t_c_br.ra.value, cord_t_c_br.dec.value])
    #         coordinates.append([cord_t_c_bl.ra.value, cord_t_c_bl.dec.value])
    #
    #         square_list.append(coordinates)

    # coord_a = SkyCoord(ra=center_ra, dec=center_dec, unit='deg')
    # # 获取方位角和高度角
    # az_alt = coord_a.transform_to(AltAz(obstime=current_time, location=location))
    # print('=========   az  [%s]    alt [%s]   ========' % (az_alt.az, az_alt.alt))
    # x_cart = coord_a.cartesian
    # print('cart  [%s]' % x_cart)
    # x_cart_array = np.array([x_cart.x.value, x_cart.y.value, x_cart.z.value])
    # left_axis, right_axis, top_axis, bottom_axis = get_l_r_t_b_axis(x_cart_array, current_time, location)
    # print('%s  %s  %s  %s  ' % (left_axis, right_axis, top_axis, bottom_axis))
    # theta = np.pi / 18
    # x_cart_l = rotate(x_cart_array, theta, left_axis)
    # x_cart_r = rotate(x_cart_array, theta, right_axis)
    # x_cart_t = rotate(x_cart_array, theta, top_axis)
    # x_cart_b = rotate(x_cart_array, theta, bottom_axis)
    # print('l:[%s]     r:[%s]' % (x_cart_l, x_cart_r))
    # print('t:[%s]     b:[%s]' % (x_cart_t, x_cart_b))
    # print('-:[%s]' % x_cart_array)
    #
    # # # 创建SkyCoord对象
    # # cartesian_coord = SkyCoord(x=x * u.pc, y=y * u.pc, z=z * u.pc, frame='gcrs', representation_type='cartesian')
    # #
    # # # 将笛卡尔坐标转换为赤经赤纬坐标
    # # icrs_coord = cartesian_coord.transform_to('gcrs')
    #
    # cart_coord_l = SkyCoord(x=x_cart_l[0] * u.pc, y=x_cart_l[1] * u.pc, z=x_cart_l[2] * u.pc,
    #                         representation_type='cartesian')
    # cart_coord_r = SkyCoord(x=x_cart_r[0] * u.pc, y=x_cart_r[1] * u.pc, z=x_cart_r[2] * u.pc,
    #                         representation_type='cartesian')
    # cart_coord_t = SkyCoord(x=x_cart_t[0] * u.pc, y=x_cart_t[1] * u.pc, z=x_cart_t[2] * u.pc,
    #                         representation_type='cartesian')
    # cart_coord_b = SkyCoord(x=x_cart_b[0] * u.pc, y=x_cart_b[1] * u.pc, z=x_cart_b[2] * u.pc,
    #                         representation_type='cartesian')
    # eq_coord_l = cart_coord_l.transform_to('gcrs')
    # eq_coord_r = cart_coord_r.transform_to('gcrs')
    # eq_coord_t = cart_coord_t.transform_to('gcrs')
    # eq_coord_b = cart_coord_b.transform_to('gcrs')
    # print('========= l:  ra  [%s]    dec [%s]   ========' % (eq_coord_l.ra, eq_coord_l.dec))
    # print('========= r:  ra  [%s]    dec [%s]   ========' % (eq_coord_r.ra, eq_coord_r.dec))
    # print('========= t:  ra  [%s]    dec [%s]   ========' % (eq_coord_t.ra, eq_coord_t.dec))
    # print('========= b:  ra  [%s]    dec [%s]   ========' % (eq_coord_b.ra, eq_coord_b.dec))
    # coordinates.append([eq_coord_l.ra.value, eq_coord_l.dec.value])
    # coordinates.append([eq_coord_r.ra.value, eq_coord_r.dec.value])
    # coordinates.append([eq_coord_t.ra.value, eq_coord_t.dec.value])
    # coordinates.append([eq_coord_b.ra.value, eq_coord_b.dec.value])
    # az_alt_l = eq_coord_l.transform_to(AltAz(obstime=current_time, location=location))
    # az_alt_r = eq_coord_r.transform_to(AltAz(obstime=current_time, location=location))
    # az_alt_t = eq_coord_t.transform_to(AltAz(obstime=current_time, location=location))
    # az_alt_b = eq_coord_b.transform_to(AltAz(obstime=current_time, location=location))
    # print('========= l:  az  [%s]    alt [%s]   ========' % (az_alt_l.az, az_alt_l.alt))
    # print('========= r:  az  [%s]    alt [%s]   ========' % (az_alt_r.az, az_alt_r.alt))
    # print('========= t:  az  [%s]    alt [%s]   ========' % (az_alt_t.az, az_alt_t.alt))
    # print('========= b:  az  [%s]    alt [%s]   ========' % (az_alt_b.az, az_alt_b.alt))
    # print('========= -:  az  [%s]    alt [%s]   ========' % (az_alt.az, az_alt.alt))
    template_root_path = "e:/test"
    output_root_path = "e:/test"
    now = datetime.now()
    time_str = current_time_with_tz.strftime('%Y%m%d_%H%M%S_%Z')
    out_path = os.path.join(output_root_path, "%s_%s.txt" % ("auto", time_str))
    # out_path = os.path.join(output_root_path, time_str, "%s_%s.txt" % ("auto", time_str))

    temp_item_file_name = "temp_item.txt"
    temp_list_file_name = "temp_list.txt"
    temp_item_file_path = os.path.join(template_root_path, temp_item_file_name)
    temp_list_file_path = os.path.join(template_root_path, temp_list_file_name)
    print("单目标模板：[%s]   计划文件模板：[%s]   输出文件:[%s]" % (temp_item_file_path, temp_list_file_path, out_path))
    print(temp_list_file_path)
    template_item_content = load_item_template(temp_item_file_path)
    write_plan_file(temp_list_file_path, center_list, out_path, template_item_content)

    return JsonResponse({'sun_ra': str(sun_radec.ra.deg), 'sun_dec': str(sun_radec.dec.deg),
                         'moon_ra': str(moon_radec.ra.deg), 'moon_dec': str(moon_radec.dec.deg),
                         'ex_message': ex_message, 'sun_curve': sun_curve, 'moon_curve': moon_curve,
                         'center_ra': str(center_ra.value), 'center_dec': str(center_dec.value),
                         'areas': square_list, 'centers': center_list, 'out_path': out_path})


def draw_load_plan(request):
    now_date = datetime.today().strftime('%Y%m%d')
    req_overlap = float(request.GET.get('req_overlap_deg', 0.2))
    req_row = int(request.GET.get('req_row_int', 2))
    req_col = int(request.GET.get('req_col_int', 2))
    req_x_img = float(request.GET.get('req_x_img_deg', 2.1))
    req_y_img = float(request.GET.get('req_y_img_deg', 4.2))
    req_ra_cen_deg = float(request.GET.get('req_ra_deg', 30.1))
    req_dec_cen_deg = float(request.GET.get('req_dec_deg', 40.3))
    ex_message = ''

    # 起始中心坐标（示例坐标，你可以根据实际需要修改）
    center_ra = req_ra_cen_deg * u.deg
    center_dec = req_dec_cen_deg * u.deg

    square_list = []
    center_list = []
    kats_region_file = "E:/test/SkyRegion.txt"
    east_region_file = "E:/test/SkyRegion_EAST.txt"
    # 打开文件
    with open(east_region_file, 'r') as file:
        # 初始化RA和DEC列表
        ra_list = []
        dec_list = []
        row_limit = 1000
        row_now = 0
        # 逐行读取文件内容
        for line in file:
            # 忽略以分号或井号开头的行
            if line.startswith(';') or line.startswith('#') or (not line.startswith('P')):
                continue

            # 提取RA和DEC
            try:
                parts = line.split()
                ra = parts[1]
                dec = parts[2]
                # 将RA和DEC添加到列表中
                ra_list.append(ra)
                dec_list.append(dec)
                ra_deg = float(ra) * 15
                center_list.append([ra_deg, dec])
                row_now = row_now + 1
                if row_now > row_limit:
                    break
                if row_now % 100 == 0:
                    print("---[%s /%s]" % (row_now, row_limit))
            except Exception as e:
                print(e)

    # # 打印提取的RA和DEC
    # for ra, dec in zip(ra_list, dec_list):
    #     print("RA:", ra, "DEC:", dec)

    num_row = req_row
    num_colum = req_col
    img_wid = req_x_img
    img_hei = req_y_img
    img_overlap = req_overlap
    row_head_coord_list = []
    cen_ra_row = center_ra
    cen_dec_row = center_dec
    coord_debug = SkyCoord(ra=center_ra, dec=center_dec, unit='deg')

    # =======================================================================================
    rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(cen_ra_row, cen_dec_row)
    rtt_l = np.array([0, 0, 1])
    rtt_r = np.array([0, 0, -1])
    # rtt_t = np.array([1, 0, 0])
    # rtt_b = np.array([-1, 0, 0])
    if num_row % 2 == 0:
        half_row = num_row // 2
        for i in range(half_row):
            cord_row_head_center_item = get_top_fix_axis(cen_ra_row, cen_dec_row, ((img_hei - img_overlap) * i)
                                                         + ((img_hei / 2) - (img_overlap / 2)), rtt_t)
            row_head_coord_list.append(cord_row_head_center_item)
        # reverse top list
        row_head_coord_list.reverse()
        for i in range(half_row):
            cord_row_head_center_item = get_bottom_fix_axis(cen_ra_row, cen_dec_row, ((img_hei - img_overlap) * i)
                                                            + ((img_hei / 2) - (img_overlap / 2)), rtt_b)
            row_head_coord_list.append(cord_row_head_center_item)
    else:
        half_row = num_row // 2
        for i in range(half_row):
            cord_row_head_center_item = get_top_fix_axis(cen_ra_row, cen_dec_row, (img_hei - img_overlap) * (i + 1),
                                                         rtt_t)
            row_head_coord_list.append(cord_row_head_center_item)
        # reverse top list
        row_head_coord_list.reverse()
        row_head_coord_list.append(SkyCoord(ra=cen_ra_row, dec=cen_dec_row, unit='deg').transform_to('gcrs'))
        for i in range(half_row):
            cord_row_head_center_item = get_bottom_fix_axis(cen_ra_row, cen_dec_row, (img_hei - img_overlap) * (i + 1),
                                                            rtt_b)
            row_head_coord_list.append(cord_row_head_center_item)

    # 计算四角点
    for i in range(len(center_list)):
        cord_real_center = center_list[i]
        rtt_l, rtt_r, rtt_t, rtt_b = get_rotate_fix_axis(cord_real_center[0], cord_real_center[1])
        # rtt_l = np.array([0, 0, 1])
        # rtt_r = np.array([0, 0, -1])
        # print('[%s]  [%s]              [%s]   [%s]' % (i, j, center_ra_colum, center_dec_colum))
        coordinates = []
        cord_t_c_tm = get_top_fix_axis(cord_real_center[0], cord_real_center[1], (img_hei / 2), rtt_t)
        cord_t_c_bm = get_bottom_fix_axis(cord_real_center[0], cord_real_center[1], (img_hei / 2), rtt_b)
        cord_t_c_tl = get_left_fix_axis(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, (img_wid / 2), rtt_l)
        cord_t_c_tr = get_right_fix_axis(cord_t_c_tm.ra.value, cord_t_c_tm.dec.value, (img_wid / 2), rtt_r)
        cord_t_c_bl = get_left_fix_axis(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, (img_wid / 2), rtt_l)
        cord_t_c_br = get_right_fix_axis(cord_t_c_bm.ra.value, cord_t_c_bm.dec.value, (img_wid / 2), rtt_r)
        coordinates.append([cord_t_c_tl.ra.value, cord_t_c_tl.dec.value])
        coordinates.append([cord_t_c_tr.ra.value, cord_t_c_tr.dec.value])
        coordinates.append([cord_t_c_br.ra.value, cord_t_c_br.dec.value])
        coordinates.append([cord_t_c_bl.ra.value, cord_t_c_bl.dec.value])

        square_list.append(coordinates)

    return JsonResponse({'ex_message': ex_message,
                         'center_ra': str(center_ra.value), 'center_dec': str(center_dec.value),
                         'areas': square_list, 'centers': center_list})


def draw_plan(request):
    plt.style.use(astropy_mpl_style)
    quantity_support()
    bj_cp = EarthLocation(lat=40 * u.deg, lon=116 * u.deg, height=390 * u.m)
    current_time = Time.now()
    m31 = SkyCoord.from_name('M31')
    utcoffset = 8 * u.hour
    midnight = Time('2023-9-23 00:00:00') - utcoffset
    delta_midnight = np.linspace(-2, 10, 100) * u.hour
    frame_July13night = AltAz(obstime=midnight + delta_midnight,
                              location=bj_cp)
    m31altazs_July13night = m31.transform_to(frame_July13night)
    delta_midnight = np.linspace(-12, 12, 100) * u.hour
    times_July12_to_13 = midnight + delta_midnight
    frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=bj_cp)
    sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)
    moon_July12_to_13 = get_body("moon", times_July12_to_13)
    moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)

    plt.plot(delta_midnight, sunaltazs_July12_to_13.alt, color='r', label='Sun')
    plt.plot(delta_midnight, moonaltazs_July12_to_13.alt, color=[0.75] * 3, ls='--', label='Moon')
    plt.scatter(delta_midnight, m31altazs_July13night.alt,
                c=m31altazs_July13night.az, label='M31', lw=0, s=8,
                cmap='viridis')
    plt.fill_between(delta_midnight, 0 * u.deg, 90 * u.deg,
                     sunaltazs_July12_to_13.alt < -0 * u.deg, color='0.5', zorder=0)
    plt.fill_between(delta_midnight, 0 * u.deg, 90 * u.deg,
                     sunaltazs_July12_to_13.alt < -18 * u.deg, color='k', zorder=0)
    plt.colorbar().set_label('Azimuth [deg]')
    plt.legend(loc='upper left')
    plt.xlim(-12 * u.hour, 12 * u.hour)
    plt.xticks((np.arange(13) * 2 - 12) * u.hour)
    plt.ylim(0 * u.deg, 90 * u.deg)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Altitude [deg]')
    # plt.show()
    # 将图像保存到内存中
    buffer = io.BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    plt.clf()
    # 返回图像作为HTTP响应
    response = HttpResponse(buffer.read(), content_type='image/png')
    buffer.close()

    return response
