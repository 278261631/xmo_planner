from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='aladin'),
    path('draw_sun', views.draw_sun, name='draw_sun'),
    path('draw_plan', views.draw_plan, name='draw_plan'),
]
