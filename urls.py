from django.conf.urls import patterns

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    # Example:
    # (r'^amview/', include('amview.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs'
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # (r'^admin/', include(admin.site.urls)),
    (r"^favicon.ico", "page_not_found"),
    (r"^(?P<path>.+/)$", "viewer.views.index"),
    (r"^(?P<path>.+/).*\.fasta$", "viewer.views.fasta"),
)
