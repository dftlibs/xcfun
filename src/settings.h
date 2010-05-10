#ifndef SETTINGS_H
#define SETTINGS_H
#include "array.h"

struct setting_data
{
  int id;
  const char *name;
  const char *short_description;  
  const char *long_description;
  double default_value;
};

class setting
{
public:
  setting() {}
  friend class user_settings;
  friend class settings_database;
protected:
  setting(const setting_data *d) : data(d) 
  { 
    id = d->id; 
  }
  const setting_data *data;
  int id;
};

class settings_database;

class user_settings
{
public:
  ~user_settings();
  friend class settings_database;
  // Set name = value, return 0 if name is a valid setting
  int set(const char *name, double value);
  double get(const char *name) const;
  double get(int n) const; // Fast lookup, unsafe
  double get(const setting &s) const; // Fast lookup, safe
  setting lookup(const char *name) const; // For use with get(setting &)
  //Functions for enumeration of possible settings
  int nr_settings(void) const;
  bool is_set(const char *name) const;
  bool is_set(int n) const ;
  const char *setting_name(int n) const;
  const char *setting_short_description(int n) const;
  const char *setting_long_description(int n) const;
  // Return the index of setting with name, or -1 if not found
  int index_of(const char *name) const;
protected:
  user_settings(const settings_database &src);
  const settings_database &db;
  double *values;
  bool *val_is_set;
  int len;
};

class settings_database
{
public:
  friend class user_settings;
  settings_database();
  ~settings_database();
  user_settings *new_user_settings(void);
  setting new_setting(const char *name, double default_value,
		      const char *short_description,
		      const char *long_description);
protected:
  bool is_frozen;
  array<const setting_data *> data;
};

#endif
