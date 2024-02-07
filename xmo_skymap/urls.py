from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='aladin'),
    path('draw_sun', views.draw_sun, name='draw_sun'),
    path('draw_all_sky_by_jgg', views.draw_all_sky_by_jgg, name='draw_all_sky_by_jgg'),
    path('draw_load_plan', views.draw_load_plan, name='draw_load_plan'),
    path('draw_plan', views.draw_plan, name='draw_plan'),
]
