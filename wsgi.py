import os
import django
import sys

sys.path.append('/srv/app/src/')
sys.path.append('/srv/app/src/amview/')

os.environ["DJANGO_SETTINGS_MODULE"] = "amview.settings"

# This application object is used by the development server
# as well as any WSGI server configured to use this file.

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()

#import django.core.handlers.wsgi
#application = django.core.handlers.wsgi.WSGIHandler()
