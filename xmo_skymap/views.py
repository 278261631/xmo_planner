import io

import numpy as np
import pytz

from astropy.coordinates import EarthLocation, AltAz, get_sun, get_moon, Angle, get_body, SkyCoord
from astropy.time import Time
from astropy.visualization import astropy_mpl_style, quantity_support
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render
from datetime import datetime
from django.core.cache import cache
from django.template import loader
import astropy.units as u
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
    ex_message = ''

    # 设置观测地点（以纽约为例）
    location = EarthLocation(lat='39.7128', lon='116.0060')  # 纬度和经度

    # 获取当前时间
    current_time = Time.now()
    current_time_with_tz = current_time.to_datetime(pytz.timezone('Asia/Shanghai'))
    # 将观测地点和时间结合起来
    altaz = AltAz(location=location, obstime=current_time)

    # 计算太阳的赤经和赤纬
    sun = get_body('sun', current_time, location)
    print(sun)
    sun_radec = sun.transform_to(AltAz(location=location, obstime=current_time)).transform_to('gcrs')

    # 计算月亮的赤经和赤纬
    moon = get_body('moon', current_time, location)
    moon_radec = moon.transform_to(AltAz(location=location, obstime=current_time)).transform_to('gcrs')
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
    print(f"当前时间：{current_time_with_tz.strftime('%Y-%m-%d %H:%M:%S %Z')}")

    return JsonResponse({'sun_ra': str(sun_radec.ra.deg), 'sun_dec': str(sun_radec.dec.deg),
                         'moon_ra': str(moon_radec.ra.deg), 'moon_dec': str(moon_radec.dec.deg),
                         'ex_message': ex_message})


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

