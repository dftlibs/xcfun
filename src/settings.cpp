#include "settings.h"
#include <cstring>

settings_database::settings_database(void)
{
  is_frozen = false;
}

settings_database::~settings_database(void)
{
  for (int i=0;i<data.size();i++)
    delete[] data[i];
}

user_settings::user_settings(const settings_database &src) : db(src)
{
  values = new double[src.data.size()];
  val_is_set = new bool[src.data.size()];
  len = src.data.size();
  for (int i=0;i<len;i++)
    {
      values[i] = src.data[i]->default_value;
      val_is_set[i] = false;
    }
}

user_settings::~user_settings()
{
  delete[] values;
  delete[] val_is_set;
}

void xc_die(const char *message, int mode);

setting user_settings::lookup(const char *name) const
{
  for (int i=0;i<db.data.size();i++)
    if (strcmp(db.data[i]->name,name) == 0)
      return setting(db.data[i]);
  xc_die("No such setting",-1);
  return setting(db.data[0]);
}


double user_settings::get(const setting &s) const
{
  return values[s.id];
}

double user_settings::get(int n) const
{
  if (n < 0 or n > len)
    return 0;
  return values[n];
}

double user_settings::get(const char *name) const
{
  int i = index_of(name);
  if (i>=0)
    return values[db.data[i]->id];
  return 0;
}

int user_settings::set(const char *name, double value)
{
  int i = index_of(name);
  if (i>=0)
    {
      values[db.data[i]->id] = value;
      val_is_set[i] = true;
      return 0;
    }
  else
    {
      return -1;
    }
}

int user_settings::nr_settings(void) const
{
  return db.data.size();
}

const char *user_settings::setting_name(int n) const
{
  assert(n>=0 and n < nr_settings());
  return db.data[n]->name;
}

const char *user_settings::setting_short_description(int n) const
{
  assert (n>=0 and n < nr_settings());
  return db.data[n]->short_description;
}

const char *user_settings::setting_long_description(int n) const
{
  assert (n>=0 and n < nr_settings());
  return db.data[n]->long_description;
}

int user_settings::index_of(const char *name) const
{
  for (int i=0;i<db.data.size();i++)
    if (strcmp(name,db.data[i]->name) == 0)
      return i;
  return -1;
}

bool user_settings::is_set(const char *name) const
{
  int i = index_of(name);
  if (i>=0)
    return val_is_set[i];
  else
    return false;
}

bool user_settings::is_set(int n) const
{
  assert (n>=0 and n < nr_settings());
  return val_is_set[n];
}

setting 
settings_database::new_setting(const char *name, double default_value,
			       const char *short_description,
			       const char *long_description)
{
  assert(not is_frozen);
  for (int i=0;i<data.size();i++)
    if (strcmp(data[i]->name,name) == 0)
      return setting(data[i]);
  setting_data *d = new setting_data;
  d->name = name;
  d->default_value = default_value;
  d->short_description = short_description;
  d->long_description = long_description;
  d->id = data.size();
  data.push_back(d);
  return setting(d);
}

user_settings *settings_database::new_user_settings(void)
{
  is_frozen = true;
  return new user_settings(*this);
}
